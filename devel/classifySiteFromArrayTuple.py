#cdarby@jhu.edu
#updated 02-01-18

#Python 2.7.13 (tested 01-24-18)
#Python 3.6.2 (tested 01-24-18)

from __future__ import print_function

import pysam,argparse,pymp,os,sys

if sys.version_info > (3,0):
    import _pickle as cPickle
else:
    import cPickle

import numpy as np
from collections import namedtuple
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib

def readIntervalFile(fname):
    R = []
    with open(fname) as bedfile:
        region = bedfile.readline().strip().split()
        while len(region) > 0:
            R.append(region)
            region = bedfile.readline().strip().split()
    return R

def readFeatures(read,featureDictShared):
    mapq = read.mapping_quality #could be None?
    if read.is_reverse: strand = 0
    else: strand = 1
    clip = read.query_length - len(read.query_alignment_sequence)
    ind = 0
    for (op,ln) in read.cigartuples:
        if op == 1 or op == 2 or op == 3:
            ind += ln

    if read.has_tag("HP"):
        XH = read.get_tag("HP")
    else:
        XH = None
    if read.has_tag("MI"): 
        MI = read.get_tag("MI")
        try:
            M = featureDictShared[MI]   
            if M[0] < 1: #initialized to 0
                XR = M[0] #nreads
                XB = M[1] #span
                XL = M[2] #read agreement fraction
                XV = M[3] #nvariants
            else:
                XR = 0
                XB = 0
                XL = 0 #or None? Change in calcLinkedReadFeatures
                XV = 0
        except IndexError:
            XR = 0
            XB = 0
            XL = 0 #or None? Change in calcLinkedReadFeatures
            XV = 0
    else: 
        XR = 0
        XB = 0
        XL = 0 #or None? Change in calcLinkedReadFeatures
        XV = 0
        
    return (mapq,strand,clip,ind,XH,XL,XV,XR,XB) 

def siteFeatures(S):
    (bases, bqs, mapqs, strands, readpos, clip, ind, XH, XL, XV, XR, XB) = (S.bases, S.bqs, S.mapqs, S.strands, S.readpos, S.clip, S.ind, S.XH, S.XL, S.XV, S.XR, S.XB)

    depth = len(bases)
    if depth == 0: 
        return None

    #features for whole site
    
    nh0 = sum([k == 0 for k in XH])
    nh1 = sum([k == 1 for k in XH])
    nphased = nh0 + nh1
    if nphased == 0: return None
    fracphased = 1.0*nphased/depth
    frach0 = 1.0*nh0/nphased
    frach1 = 1.0*nh1/nphased
    frach = max(frach0,frach1)
    fracstrand = 1.0*sum(strands)/depth
    
    values, counts = np.unique(bases,return_counts=True)
    if len(values) == 1: return None #no minor allele bases
    basecounts = zip(values,counts)
    bases_ranked = sorted(basecounts, key=lambda x:x[1])
    major = bases_ranked[-1][0] #most frequent nt
    minor = bases_ranked[-2][0] #if bases_ranked[-2][0] != "N" else bases_ranked[-3][0]
    nmajor = bases_ranked[-1][1]
    nminor = bases_ranked[-2][1] #if bases_ranked[-2][0] != "N" else bases_ranked[-3][1]
    MAF = 1.0*nminor/depth
    fracOtherAllele = 1.0*(depth - nmajor - nminor)/depth
    MAFnorm = abs(0.25-MAF)
    MAF_phased = 0.0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if bases[i] != major and XH[i] != None])/nphased #MAF of phased reads

    #linked features
    #Reads that need to be changed in the allele x haplotype matrix
    change_to_hom_mask = [True if (XH[i] != None and bases[i] == minor) else False for i in range(depth)] #ignore non major/minor bases
    change_to_het_left_mask = [True if (XH[i] == 1 and bases[i] == major) or (XH[i] == 0 and bases[i] == minor) else False for i in range(depth) ]
    change_to_het_right_mask = [True if (XH[i] == 0 and bases[i] == major) or (XH[i] == 1 and bases[i] == minor) else False for i in range(depth) ]

    change_to_hom = sum(change_to_hom_mask)
    change_to_het_right = sum(change_to_het_right_mask)
    change_to_het_left = sum(change_to_het_left_mask)
    if change_to_hom <= change_to_het_left and change_to_hom <= change_to_het_right:
        mask = change_to_hom_mask
    elif change_to_het_left <= change_to_hom  and change_to_het_left <= change_to_het_right:
        mask = change_to_het_left_mask
    else:   
        mask = change_to_het_right_mask

    nC = sum(mask)
    if nC == 0 or nC == nphased: #No C-reads or no N-reads
        return None 

    Cavgbq = 1.0*np.compress(mask,bqs).sum()/nC
    C_haps = np.compress(mask,XH)
    Cfrach0 = 1.0*sum([1 for i in range(nC) if C_haps[i] == 0])/nC
    Cfrach1 = 1.0*sum([1 for i in range(nC) if C_haps[i] == 1])/nC
    Cfrach = max(Cfrach0,Cfrach1)
    
    Cavgmapq = 1.0*np.compress(mask,mapqs).sum()/nC
    Cfracstrand = 1.0*np.compress(mask,strands).sum()/nC
    C_bases = np.compress(mask,bases)
    values, counts = np.unique(C_bases,return_counts=True)
    basecounts = zip(values,counts)
    bases_ranked = sorted(basecounts, key=lambda x:x[1])
    CMAF = 1.0*bases_ranked[-1][1]/depth
    
    CavgXL = 1.0*np.compress(mask,XL).sum()/nC
    CavgXV = 1.0*np.compress(mask,XV).sum()/nC
    CavgXR = 1.0*np.compress(mask,XR).sum()/nC
    CavgXB = 1.0*np.compress(mask,XB).sum()/nC
    Cavgpos = 1.0*np.compress(mask,readpos).sum()/nC
    Cavgclip = 1.0*np.compress(mask,clip).sum()/nC
    Cavgind = 1.0*np.compress(mask,ind).sum()/nC

    mask_inv = [False if (mask[i] or XH[i] == None or bases[i] != major or bases[i] != minor) else True for i in range(depth) ]
    nN = sum(mask_inv)
    Navgbq = 1.0*np.compress(mask_inv,bqs).sum()/nN
    N_haps = np.compress(mask_inv,XH)
    Nfrach0 = 1.0*sum([1 for i in range(nN) if N_haps[i] == 0])/nN
    Nfrach1 = 1.0*sum([1 for i in range(nN) if N_haps[i] == 1])/nN
    Nfrach = max(Nfrach0,Nfrach1)
    
    Navgmapq = 1.0*np.compress(mask_inv,mapqs).sum()/nN
    Nfracstrand = 1.0*np.compress(mask_inv,strands).sum()/nN
    
    N_bases = np.compress(mask_inv,bases)
    values, counts = np.unique(N_bases,return_counts=True)
    basecounts = zip(values,counts)
    bases_ranked = sorted(basecounts, key=lambda x:x[1])
    NMAF = 1.0*bases_ranked[-1][1]/depth

    NavgXL = 1.0*np.compress(mask,XL).sum()/nN
    NavgXV = 1.0*np.compress(mask,XV).sum()/nN
    NavgXR = 1.0*np.compress(mask,XR).sum()/nN
    NavgXB = 1.0*np.compress(mask,XB).sum()/nN
    Navgpos = 1.0*np.compress(mask,readpos).sum()/nN
    Navgclip = 1.0*np.compress(mask,clip).sum()/nN
    Navgind = 1.0*np.compress(mask,ind).sum()/nN

    fracC = 1.0*nC/(nN+nC)
    weightedCbq = 1.0*nC*Cavgbq/(nC*Cavgbq+nN*Navgbq)
    Cavgibq = 1-weightedCbq #Should rm this feature - redundant. #1.0*sum(1-x for x in Nbqs)/sum([1-x for x in Cbqs+Nbqs])
    
    #non-linked features - "M-reads" are those with minor allele (like C-reads), "J-reads" are those with major allele (like N-reads)
    mask_M = [True if bases[i] == minor else False for i in range(depth)]
    Mavgbq = 1.0*np.compress(mask_M,bqs).sum()/nminor
    Mavgmapq = 1.0*np.compress(mask_M,mapqs).sum()/nminor
    Mfracstrand = 1.0*np.compress(mask_M,strands).sum()/nminor
    Mavgpos = 1.0*np.compress(mask_M,readpos).sum()/nminor
    Mavgclip = 1.0*np.compress(mask_M,clip).sum()/nminor
    Mavgind = 1.0*np.compress(mask_M,ind).sum()/nminor
    
    mask_J = [True if bases[i] == major else False for i in range(depth)]
    Javgbq = 1.0*np.compress(mask_J,bqs).sum()/nmajor
    Javgmapq = 1.0*np.compress(mask_J,mapqs).sum()/nmajor
    Jfracstrand = 1.0*np.compress(mask_J,strands).sum()/nmajor
    Javgpos = 1.0*np.compress(mask_J,readpos).sum()/nmajor
    Javgclip = 1.0*np.compress(mask_J,clip).sum()/nmajor
    Javgind = 1.0*np.compress(mask_J,ind).sum()/nmajor
    
    weightedMbq = 1.0*np.compress(mask_M,bqs).sum()/(Mavgbq * nminor + Javgbq * nmajor)
    Mavgibq = 1-weightedMbq #Should rm this feature - redundant. 
    features = np.array([depth,fracphased,fracstrand,frach,MAF,MAF_phased,
        fracC,Cavgbq,Cfrach,Cavgmapq,Cfracstrand,CMAF,CavgXL,CavgXV,CavgXR,CavgXB,Cavgpos, Cavgclip,Cavgind,
        Navgbq,Nfrach,Navgmapq,Nfracstrand,NMAF,NavgXL,NavgXV,NavgXR,NavgXB,Navgpos,Navgclip,Navgind,
        weightedCbq,Cavgibq,
        Mavgbq,Mavgmapq,Mfracstrand,Mavgpos,Mavgclip,Mavgind,Mavgibq,weightedMbq,
        Javgbq,Javgmapq,Jfracstrand,Javgpos,Javgclip,Javgind,MAFnorm,fracOtherAllele])
    return features



parser = argparse.ArgumentParser(description='python classifySiteFromDictTuple.py')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name',required=True)
parser.add_argument('--featuredict', help='.pkl file from calcLinkedReadFeatures.py' ,required=True)
parser.add_argument('--clf',help="random forest scikit-learn classifier",required=True)
parser.add_argument('--nproc',help="parallelism",required=False,default=3)

args = parser.parse_args()
R = readIntervalFile(args.bed)
Site = namedtuple('Site','pos bases bqs mapqs strands readpos clip ind XH XL XV XR XB')


NPAR = int(args.nproc)
READBATCH = 100000
CLFBATCH = 10000
PRINT_TH = 0.8
NFEAT = 49

featureDict = cPickle.load(open(args.featuredict,'rb'))
NROWS = max(featureDict.keys()) + 1
featureArr = pymp.shared.array((NROWS,4),dtype='float32') #default float64
for (k,v) in featureDict.items():
    featureArr[k][0] = v[0]
    featureArr[k][1] = v[1]
    featureArr[k][2] = v[2]
    featureArr[k][3] = v[3]
featureDict = None #so that it's not copied by each parallel

clf = joblib.load(args.clf)
samfile = pysam.AlignmentFile(args.bam, "rb") #could do this per-parallel if bamfile is split


with pymp.Parallel(NPAR) as pymp_context:
    for region_index in pymp_context.xrange(len(R)): #pymp.xrange returns an iterator and corresponds to dynamic scheduling.
        #pymp_context.print(pymp_context.thread_num)
        region = R[region_index]     
        chrom = region[0]
        start = int(region[1])
        end = int(region[2])
        predictionArray = np.empty([CLFBATCH,NFEAT]) #batch of feature vectors ready for prediction
        positionsPredicted = [] #bp positions of rows in predictionArray
        nPred = 0
        SITES = dict()
        last_site_eval = -1
        nreads = 0

        for read in samfile.fetch(chrom, start, end, multiple_iterators=True): 
            #Each iteration returns a AlignedSegment object which represents a single read along with its fields and optional tags - Returns in order by start alignment position
            if read.is_duplicate or read.is_qcfail: continue #or aln.is_secondary or aln.is_unmapped or aln.is_supplementary
            #reject secondary or supplementary?
            nreads += 1
            
            #Features for the whole read
            (mapq,strand,clip,ind,XH,XL,XV,XR,XB) = readFeatures(read,featureArr)#Shared)

            aligned_pairs = read.get_aligned_pairs()
            for p in aligned_pairs:
                if p[0] == None or p[1] == None or read.query_sequence[p[0]] == "N":
                    continue
                querypos = p[0]
                refpos = p[1]

                if refpos in SITES:
                    S = SITES[refpos]
                else:
                    #S = Site(refpos)
                    S = Site(refpos,[],[],[],[],[],[],[],[],[],[],[],[])
                    SITES[refpos] = S

                S.bases.append(read.query_sequence[querypos])
                S.bqs.append(read.query_qualities[querypos])
                S.mapqs.append(mapq)
                S.strands.append(strand)
                S.readpos.append(min((read.query_alignment_end - querypos),(querypos - read.query_alignment_start)))
                S.clip.append(clip)
                S.ind.append(ind)
                S.XH.append(XH)
                S.XL.append(XL)
                S.XV.append(XV)
                S.XR.append(XR)
                S.XB.append(XB)

            if nreads >= READBATCH: #batch size for evaluating sites where all reads seen 
                #pymp_context.print("Read batch")
                #pymp_context.print(pymp_context.thread_num)
                #pymp_context.print(sys.getsizeof(SITES))
                #pymp_context.print(min(read.reference_start,end)-max(last_site_eval+1,start))
                L = len(SITES)
                for refpos in range(max(last_site_eval+1,start), min(read.reference_start,end)):
                    if refpos in SITES:
                        S = SITES[refpos]
                    else: continue #eek!

                    feature_vector = siteFeatures(S)
                    if type(feature_vector) != type(None):
                        predictionArray[nPred] = feature_vector
                        nPred += 1
                        positionsPredicted.append(refpos)
                    del SITES[S.pos]

                    if nPred == CLFBATCH:
                        #pymp_context.print("Decision tree")
                        #pymp_context.print(sys.getsizeof(predictionArray))
                        p = clf.predict_proba(predictionArray)
                        for i in range(len(p)):
                            if p[i][1] >= PRINT_TH:
                                #continue
                                pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + " " + str(round(p[i][1],3))+ " " + str(round(predictionArray[i][4],3)) + " " + str(predictionArray[i][0]) + " " + str(round(predictionArray[i][1],3)))
                        positionsPredicted = [] #bp positions of rows in predictionArray
                        nPred = 0
                #pymp_context.print("After read batch")
                #pymp_context.print(L-len(SITES))
                last_site_eval = read.reference_start - 1
                nreads = 0
        #pymp_context.print(sys.getsizeof(SITES))
        for (refpos,S) in SITES.items(): #After batches, finish remaining reads
            if refpos < start or refpos >= end: continue
            feature_vector = siteFeatures(S)
            if type(feature_vector) != type(None):
                predictionArray[nPred] = feature_vector
                nPred += 1
                positionsPredicted.append(refpos)

            if nPred == CLFBATCH:
                #pymp_context.print("Nearly last predict")
                p = clf.predict_proba(predictionArray)
                for i in range(len(p)):
                    if p[i][1] >= PRINT_TH:
                        #continue
                        pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + " " + str(round(p[i][1],3))+ " " + str(round(predictionArray[i][4],3)) + " " + str(predictionArray[i][0]) + " " + str(round(predictionArray[i][1],3)))
                positionsPredicted = [] #bp positions of rows in predictionArray
                nPred = 0
        #pymp_context.print("Final predict")
        #pymp_context.print(nPred)
        if nPred > 0: #After last partial-batch of reads, clf remaining feature vectors
            p = clf.predict_proba(predictionArray[:nPred])
            for i in range(nPred):
                if p[i][1] >= PRINT_TH:
                        #continue
                        pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + " " + str(round(p[i][1],3))+ " " + str(round(predictionArray[i][4],3)) + " " + str(predictionArray[i][0]) + " " + str(round(predictionArray[i][1],3)))
        pymp_context.print("*" + chrom + "*" + str(start) + "*" + str(end) + "*")

