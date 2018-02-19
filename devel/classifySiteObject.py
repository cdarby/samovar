#cdarby@jhu.edu
#updated 01-22-18

from __future__ import print_function

import pysam,argparse,pymp,os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib


NT = "ACGTN"
NPAR = 3
READBATCH = 1000000
CLFBATCH = 100000
PRINT_TH = 0
NFEAT = 49

def readIntervalFile(fname):
    R = []
    with open(fname) as bedfile:
        region = bedfile.readline().strip().split()
        while len(region) > 0:
            R.append(region)
            region = bedfile.readline().strip().split()
    return R


class Site:
    def __init__(self, pos):
        self.pos = pos
        self.bases = []
        self.bqs = []
        self.mapqs = []
        self.strands = []
        self.readpos= []
        self.clip = []
        self.ind = []
        self.XH = []
        self.XL = []
        self.XV = []
        self.XR = []
        self.XB = []       


def readFeatures(read):
    mapq = read.mapping_quality #could be None?
    if read.is_reverse: strand = 0
    else: strand = 1
    clip = read.query_length - len(read.query_alignment_sequence)
    ind = 0
    for (op,ln) in read.cigartuples:
        if op == 1 or op == 2 or op == 3:
            ind += ln

    if read.has_tag("XH"): XH = read.get_tag("XH")
    else: XH = None
    if read.has_tag("XL"): XL = read.get_tag("XL")
    else: XL = None
    if read.has_tag("XV"):
        XV = read.get_tag("XV")
        XR = read.get_tag("XR")
        XB = read.get_tag("XB")
    else:
        XV = 0
        XR = 0
        XB = 0
    return (mapq,strand,clip,ind,XH,XL,XV,XR,XB) 

def siteFeatures(S):
    (bases, bqs, mapqs, strands, readpos, clip, ind, XH, XL, XV, XR, XB) = (np.array(S.bases), np.array(S.bqs), np.array(S.mapqs), np.array(S.strands), np.array(S.readpos), np.array(S.clip), np.array(S.ind), np.array(S.XH), np.array(S.XL), np.array(S.XV), np.array(S.XR), np.array(S.XB))

    depth = len(bases)
    if depth == 0: 
        return None

    #FROM USUAL PROCEDURE
    #features for whole site
    
    nh0 = (XH == 0).sum()
    nh1 = (XH == 1).sum()
    nphased = nh0 + nh1
    if nphased == 0: return None
    fracphased = 1.0*nphased/depth
    frach0 = 1.0*nh0/nphased
    frach1 = 1.0*nh1/nphased
    frach = max(frach0,frach1)
    fracstrand = 1.0*strands.sum()/depth
    
    values, counts = np.unique(bases,return_counts=True)
    if len(values) == 1: return None #no minor allele bases
    basecounts = zip(values,counts)
    bases_ranked = sorted(basecounts, key=lambda x:x[1])
    major = bases_ranked[-1][0] #most frequent nt
    minor = bases_ranked[-2][0] #if bases_ranked[-2][0] != "N" else bases_ranked[-3][0]
    nmajor = bases_ranked[-1][1]
    nminor = bases_ranked[-2][1] #if bases_ranked[-2][0] != "N" else bases_ranked[-3][1]
    MAF = 1.0*bases_ranked[-1][1]/depth
    fracOtherAllele = 1.0*(depth - nmajor - nminor)/depth
    MAFnorm = abs(0.25-MAF)
    MAF_phased = 0.0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if bases[i] != major and XH[i] != None])/nphased #MAF of phased reads

    #linked features
    #Reads that need to be changed in the allele x haplotype matrix
    change_to_hom_mask = [True if (XH[i] != None and bases[i] != major) else False for i in range(depth)]
    change_to_het_left_mask = [True if (XH[i] == 1 and bases[i] == major) or (XH[i] == 0 and bases[i] != major) else False for i in range(depth) ]
    change_to_het_right_mask = [True if (XH[i] == 0 and bases[i] == major) or (XH[i] == 1 and bases[i] != major) else False for i in range(depth) ]

    change_to_hom = sum(change_to_hom_mask)
    change_to_het_right = sum(change_to_het_right_mask)
    change_to_het_left = sum(change_to_het_left_mask)
    if change_to_hom <= change_to_het_left and change_to_hom <= change_to_het_right:
        mask = change_to_hom_mask
    elif change_to_het_left <= change_to_hom  and change_to_het_left <= change_to_het_right:
        mask = change_to_het_left_mask
    else:   
        mask = change_to_het_right_mask

    if sum(mask) == 0 or sum(mask) == nphased: #No C-reads or no N-reads
        return None 

    nC = sum(mask)
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

    mask_inv = [False if (mask[i] or XH[i] == None) else True for i in range(depth) ]
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
    Javgbq = 1.0*np.compress(mask_J,bqs).sum()/nminor
    Javgmapq = 1.0*np.compress(mask_J,mapqs).sum()/nminor
    Jfracstrand = 1.0*np.compress(mask_J,strands).sum()/nminor
    Javgpos = 1.0*np.compress(mask_J,readpos).sum()/nminor
    Javgclip = 1.0*np.compress(mask_J,clip).sum()/nminor
    Javgind = 1.0*np.compress(mask_J,ind).sum()/nminor
    weightedMbq = 1.0*np.compress(mask_M,bqs).sum()/(Mavgbq * nminor + Javgbq * nmajor)
    Mavgibq = 1-weightedMbq #Should rm this feature - redundant. 
    features = np.array([depth,fracphased,fracstrand,frach,MAF,MAF_phased,
        fracC,Cavgbq,Cfrach,Cavgmapq,Cfracstrand,CMAF,CavgXL,CavgXV,CavgXR,CavgXB,Cavgpos, Cavgclip,Cavgind,
        Navgbq,Nfrach,Navgmapq,Nfracstrand,NMAF,NavgXL,NavgXV,NavgXR,NavgXB,Navgpos,Navgclip,Navgind,
        weightedCbq,Cavgibq,
        Mavgbq,Mavgmapq,Mfracstrand,Mavgpos,Mavgclip,Mavgind,Mavgibq,weightedMbq,
        Javgbq,Javgmapq,Jfracstrand,Javgpos,Javgclip,Javgind,MAFnorm,fracOtherAllele])
    return features



parser = argparse.ArgumentParser(description='python siteFeatures.py --bam diploid_output/phased_possorted_bam.bam --bed haploidPos.bed --out diploidBam_haploidPos.pileup')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name',required=True)
#parser.add_argument('--out', help='Output file name',required=True)
#parser.add_argument('--offset',help="Put the offset from the bed file in case an indexing mistake occurred",required=False,default=0)
#parser.add_argument('--avoid',help="sites to avoid, format chr1\t31553290, file names split with comma",required=False,default=None)
parser.add_argument('--clf',help="random forest scikit-learn classifier",required=True)

args = parser.parse_args()

R = readIntervalFile(args.bed)

with pymp.Parallel(NPAR) as pymp_context:
    for region_index in pymp_context.range(len(R)):
        region = R[region_index]     
        chrom = region[0]
        samfile = pysam.AlignmentFile(args.bam, "rb") #each thread opens file - may depend on chrom if file is split
        clf = joblib.load(args.clf) #each thread loads clf

        start = int(region[1])
        end = int(region[2])
        predictionArray = np.empty([CLFBATCH,NFEAT]) #batch of feature vectors ready for prediction
        positionsPredicted = [] #bp positions of rows in predictionArray
        nPred = 0
        SITES = dict()
        last_site_eval = -1
        nreads = 0

        for read in samfile.fetch(chrom, start, end): 
            #Each iteration returns a AlignedSegment object which represents a single read along with its fields and optional tags - Returns in order by start alignment position
            if read.is_duplicate: continue
            #reject secondary or supplementary?
            nreads += 1
            
            #Features for the whole read
            (mapq,strand,clip,ind,XH,XL,XV,XR,XB) = readFeatures(read)

            aligned_pairs = read.get_aligned_pairs()
            for p in aligned_pairs:
                if p[0] == None or p[1] == None or read.query_sequence[p[0]] == "N":
                    continue
                querypos = p[0]
                refpos = p[1]

                if refpos in SITES:
                    S = SITES[refpos]
                else:
                    S = Site(refpos)
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
                        p = clf.predict_proba(predictionArray)
                        print(p)
                        for i in range(len(p)):
                            if p[i][1] >= PRINT_TH:
                                pass
                                #pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + " " + str(round(p[i][1],3))+ " " + str(round(predictionArray[i][4],3)) + " " + str(predictionArray[i][0]) + " " + str(round(predictionArray[i][1],3)))
                        positionsPredicted = [] #bp positions of rows in predictionArray
                        nPred = 0

                last_site_eval = read.reference_start - 1
                nreads = 0

        for (refpos,S) in SITES.items(): #After batches, finish remaining reads
            if refpos < start or refpos >= end: continue
            feature_vector = siteFeatures(S)
            if type(feature_vector) != type(None):
                predictionArray[nPred] = feature_vector
                nPred += 1
                positionsPredicted.append(refpos)

            if nPred == CLFBATCH:
                #print("Nearly last predict")
                p = clf.predict_proba(predictionArray)
                for i in range(len(p)):
                    if p[i][1] >= PRINT_TH:
                        pass
                        #pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + " " + str(round(p[i][1],3))+ " " + str(round(predictionArray[i][4],3)) + " " + str(predictionArray[i][0]) + " " + str(round(predictionArray[i][1],3)))
                positionsPredicted = [] #bp positions of rows in predictionArray
                nPred = 0
        #print("Final predict")
        #print(nPred)
        if nPred > 0: #After last partial-batch of reads, clf remaining feature vectors
            p = clf.predict_proba(predictionArray)
            for i in range(nPred):
                if p[i][1] >= PRINT_TH:
                        pass
                        #pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + " " + str(round(p[i][1],3))+ " " + str(round(predictionArray[i][4],3)) + " " + str(predictionArray[i][0]) + " " + str(round(predictionArray[i][1],3)))

