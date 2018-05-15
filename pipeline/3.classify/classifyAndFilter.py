#cdarby@jhu.edu
#updated 02-02-18

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


##### Filter Parameters #####
MIN_DEPTH = 16
MIN_FRAC_PHASED = 0.5
MIN_FRAC_HAP1 = 0.3
MAX_FRAC_OTHERALLELE = 0.05
MIN_MAF = 0.05
MIN_HAPDISCORD_READS = 4
MIN_HAPDISCORD_READS_HAP = 0.1
MIN_AVG_POS_FROM_END = 10 #not for haplotype-concordant reads
MAX_AVG_INDELS = 5
MIN_CLF_SCORE = 0.05

SOFT_CLIP = True
GOOD_PAIR = True
STRAND = True

##### Performance parameters ####
READBATCH = 50000
CLFBATCH = 1000
INTERVALBATCH = 8000

def readFeatures(read):
    clipLen = read.query_length - len(read.query_alignment_sequence)
    if clipLen > 1: clip = 1
    else: clip = 0
    ind = 0
    C = read.cigartuples
    if C is not None: #When is it None? 
        for (op,ln) in C:
            if op == 1 or op == 2 or op == 3:
                ind += ln
    ASXS = read.get_tag("AS") - read.get_tag("XS")
    if read.has_tag("HP"):
        HP = int(read.get_tag("HP"))-1
    else:
        HP = None  
    if read.is_proper_pair: pair = 0
    else: pair = 1
    if read.is_reverse: strands = 0
    else: strands = 1      
    
    return (clip,clipLen,ind,HP,ASXS,pair,strands) 

def siteFeatures(S):
    (bases, bqs, readpos, clip, ind, HP, ASXS, pair, strands, clipLen) = (S.bases, S.bqs, S.readpos, S.clip, S.ind, S.HP, S.ASXS, S.pair, S.strands, S.clipLen)

    depth = len(bases)
    if depth < MIN_DEPTH: 
        #pymp_context.print("Depth")
        return None

    #features for whole site
    
    nh0 = sum([k == 0 for k in HP])
    nh1 = sum([k == 1 for k in HP])
    nphased = nh0 + nh1
    if nphased < MIN_HAPDISCORD_READS: 
        return None
    fracphased = 1.0*nphased/depth
    if fracphased < MIN_FRAC_PHASED: 
        #pymp_context.print("Fraction phased " + str(fracphased))
        return None
    frach0 = 1.0*nh0/nphased
    frach1 = 1.0*nh1/nphased
    frach = max(frach0,frach1)
    if frach < MIN_FRAC_HAP1 or (1-frach) < MIN_FRAC_HAP1: 
        #pymp_context.print("Fraction haplotype " + str(frach))
        return None
    
    values, counts = np.unique(bases,return_counts=True)
    if len(values) == 1: return None #no minor allele bases
    basecounts = zip(values,counts)
    bases_ranked = sorted(basecounts, key=lambda x:x[1])
    major = bases_ranked[-1][0] #most frequent nt
    minor = bases_ranked[-2][0] #if bases_ranked[-2][0] != "N" else bases_ranked[-3][0]
    nmajor = bases_ranked[-1][1]
    nminor = bases_ranked[-2][1] #if bases_ranked[-2][0] != "N" else bases_ranked[-3][1]
    MAF = 1.0*nminor/depth
    MAF_phased = 0.0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if bases[i] != major and HP[i] != None])/nphased
    if MAF < MIN_MAF: return None
    fracOtherAllele = 1.0*(depth - nmajor - nminor)/depth
    if fracOtherAllele > MAX_FRAC_OTHERALLELE: return None

    #linked features
    #Reads that need to be changed in the allele x haplotype matrix
    change_to_hom_mask = [True if (HP[i] != None and bases[i] == minor) else False for i in range(depth)] #ignore non major/minor bases
    change_to_het_left_mask = [True if (HP[i] == 1 and bases[i] == major) or (HP[i] == 0 and bases[i] == minor) else False for i in range(depth) ]
    change_to_het_right_mask = [True if (HP[i] == 0 and bases[i] == major) or (HP[i] == 1 and bases[i] == minor) else False for i in range(depth) ]

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
    if nC < MIN_HAPDISCORD_READS or nC == nphased: return None 
    C_haps = np.compress(mask,HP)
    Cfrach0 = 1.0*sum([1 for i in range(nC) if C_haps[i] == 0])/nC
    Cfrach1 = 1.0*sum([1 for i in range(nC) if C_haps[i] == 1])/nC
    Cfrach = max(Cfrach0,Cfrach1)
    if not (Cfrach < MIN_HAPDISCORD_READS_HAP or (1-Cfrach) < MIN_HAPDISCORD_READS_HAP): 
        #pymp_context.print("Number C-reads filter")
        return None

    Cavgpos = 1.0*np.compress(mask,readpos).sum()/nC #average distance from end of alignments
    Cavgind = 1.0*np.compress(mask,ind).sum()/nC #average number of indels in alignment
    Cclip = np.compress(mask,clip).sum() #number of reads clipped
    Cpair = np.compress(mask,pair).sum() #number of reads with bad pair
    Cstrand = np.compress(mask,strands).sum() #number of reads on + strand
    if Cavgpos < MIN_AVG_POS_FROM_END or Cavgind > MAX_AVG_INDELS or (SOFT_CLIP and Cclip == nC) or (GOOD_PAIR and Cpair == nC) or (STRAND and (Cstrand == 0 or Cstrand == nC)): 
        #pymp_context.print("C-read filter")
        return None

    mask_inv = [False if (mask[i] or HP[i] == None or (bases[i] != major and bases[i] != minor)) else True for i in range(depth) ]
    nN = sum(mask_inv)

    Navgpos = 1.0*np.compress(mask_inv,readpos).sum()/nN
    Navgind = 1.0*np.compress(mask_inv,ind).sum()/nN
    Npair = np.compress(mask_inv,pair).sum()
    Nclip = np.compress(mask_inv,clip).sum() #number of reads clipped
    Nstrand = np.compress(mask_inv,strands).sum() #number of reads on + strand
    if Navgind > MAX_AVG_INDELS or (SOFT_CLIP and Nclip == nN) or (GOOD_PAIR and Npair == nN) or (STRAND and (Nstrand == 0 or Nstrand == nN)): 
        #pymp_context.print("N-read filter")
        return None

    #non-linked features - "M-reads" are those with minor allele (like C-reads), "J-reads" are those with major allele (like N-reads)
    mask_M = [True if bases[i] == minor else False for i in range(depth)]
    Mavgpos = 1.0*np.compress(mask_M,readpos).sum()/nminor
    Mavgind = 1.0*np.compress(mask_M,ind).sum()/nminor
    Mclip = np.compress(mask_M,clip).sum()
    Mpair = np.compress(mask_M,pair).sum()
    Mstrand = np.compress(mask_M,strands).sum()
    if Mavgpos < MIN_AVG_POS_FROM_END or Mavgind > MAX_AVG_INDELS or (SOFT_CLIP and Mclip == nminor) or (GOOD_PAIR and Mpair == nminor) or (STRAND and (Mstrand == 0 or Mstrand == nminor)): 
        #pymp_context.print("M-read filter")
        return None

    mask_J = [True if bases[i] == major else False for i in range(depth)]
    Javgpos = 1.0*np.compress(mask_J,readpos).sum()/nmajor
    Javgind = 1.0*np.compress(mask_J,ind).sum()/nmajor
    Jclip = np.compress(mask_J,clip).sum()
    Jpair = np.compress(mask_J,pair).sum()
    Jstrand = np.compress(mask_J,strands).sum()
    if Javgpos < MIN_AVG_POS_FROM_END or Javgind > MAX_AVG_INDELS or (SOFT_CLIP and Jclip == nmajor) or (GOOD_PAIR and Jpair == nmajor) or (STRAND and (Jstrand == 0 or Jstrand == nmajor)): 
        #pymp_context.print("J-read filter")
        return None

    #if passes filters, do these computations

    C_bases = np.compress(mask,bases)
    values, counts = np.unique(C_bases,return_counts=True)
    basecounts = zip(values,counts)
    bases_ranked = sorted(basecounts, key=lambda x:x[1])
    CMAF = 1.0*bases_ranked[-1][1]/depth
    CavgASXS = 1.0*np.compress(mask,ASXS).sum()/nC
    Cavgclip = 1.0*np.compress(mask,clipLen).sum()/nC #average number of bases clipped
    fracC = 1.0*nC/(nN+nC)
    Cavgbq = 1.0*np.compress(mask,bqs).sum()/nC
    

    N_bases = np.compress(mask_inv,bases)
    values, counts = np.unique(N_bases,return_counts=True)
    basecounts = zip(values,counts)
    bases_ranked = sorted(basecounts, key=lambda x:x[1])
    NMAF = 1.0*bases_ranked[-1][1]/depth
    NavgASXS = 1.0*np.compress(mask_inv,ASXS).sum()/nN
    Navgclip = 1.0*np.compress(mask_inv,clipLen).sum()/nN

    Navgbq = 1.0*np.compress(mask_inv,bqs).sum()/nN
    N_haps = np.compress(mask_inv,HP)
    Nfrach0 = 1.0*sum([1 for i in range(nN) if N_haps[i] == 0])/nN
    Nfrach1 = 1.0*sum([1 for i in range(nN) if N_haps[i] == 1])/nN
    Nfrach = max(Nfrach0,Nfrach1)
    weightedCbq = 1.0*nC*Cavgbq/(nC*Cavgbq+nN*Navgbq)

    Mavgbq = 1.0*np.compress(mask_M,bqs).sum()/nminor
    Mavgclip = 1.0*np.compress(mask_M,clipLen).sum()/nminor
    MavgASXS = 1.0*np.compress(mask_M,ASXS).sum()/nminor
    
    Javgbq = 1.0*np.compress(mask_J,bqs).sum()/nmajor
    Javgclip = 1.0*np.compress(mask_J,clipLen).sum()/nmajor
    JavgASXS = 1.0*np.compress(mask_J,ASXS).sum()/nmajor
    weightedMbq = 1.0*np.compress(mask_M,bqs).sum()/(Mavgbq * nminor + Javgbq * nmajor)

    #features for classifier
    features = np.array([depth,fracphased,frach,MAF,MAF_phased,
        nC,fracC,Cavgbq,Cfrach,CMAF,
        Cavgpos,Cavgclip,Cavgind,CavgASXS,Navgbq,
        Nfrach,NMAF,Navgpos,Navgclip,Navgind,
        NavgASXS,weightedCbq,Mavgbq,Mavgpos,Mavgclip,
        Mavgind,MavgASXS,weightedMbq,Javgbq,Javgpos,
        Javgclip,Javgind,JavgASXS])
    return features

parser = argparse.ArgumentParser(description='python classify.py --bam sample.bam --bed genome.bed --clf clf.pkl --nproc 8')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name',required=True)
parser.add_argument('--clf',help="random forest scikit-learn classifier",required=True)
parser.add_argument('--nproc',help="parallelism",required=False,default=3)

args = parser.parse_args()
Site = namedtuple('Site','pos bases bqs readpos clip clipLen ind HP ASXS pair strands')

NFEAT = 33
NPAR = int(args.nproc)

clf = joblib.load(args.clf)
samfile = pysam.AlignmentFile(args.bam, "rb") #could do this per-parallel if bamfile is split

with open(args.bed) as bedfile:
    r = bedfile.readline().strip().split()
    while len(r) > 0:
        i = 0
        R = []
        while i < INTERVALBATCH and len(r) > 0:
            i+=1
            R.append(r)
            r = bedfile.readline().strip().split()

        with pymp.Parallel(NPAR,if_=(NPAR > 1)) as pymp_context:
            for region_index in pymp_context.xrange(len(R)): #pymp.xrange returns an iterator and corresponds to dynamic scheduling.
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
                    if read.is_duplicate or read.is_qcfail or read.is_secondary or read.is_supplementary: continue # or aln.is_unmapped
                    nreads += 1
                    
                    #Features for the whole read
                    (clip,clipLen,ind,HP,ASXS,pair,strands) = readFeatures(read)

                    aligned_pairs = read.get_aligned_pairs()
                    for p in aligned_pairs:
                        if p[0] == None or p[1] == None:# or read.query_sequence[p[0]] == "N":
                            continue
                        querypos = p[0]
                        refpos = p[1]

                        if refpos in SITES:
                            S = SITES[refpos]
                        else:
                            S = Site(refpos,[],[],[],[],[],[],[],[],[],[])
                            SITES[refpos] = S

                        S.bases.append(read.query_sequence[querypos])
                        S.bqs.append(read.query_qualities[querypos])
                        S.readpos.append(min((read.query_alignment_end - querypos),(querypos - read.query_alignment_start)))
                        S.clip.append(clip)
                        S.ind.append(ind)
                        S.HP.append(HP)
                        S.ASXS.append(ASXS)
                        S.pair.append(pair)
                        S.strands.append(strands)
                        S.clipLen.append(clipLen)

                    if nreads >= READBATCH: #batch size for evaluating sites where all reads seen 
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
                                #pymp_context.print(str(feature_vector))
                            del SITES[S.pos]

                            if nPred == CLFBATCH:
                                p = clf.predict_proba(predictionArray)
                                for i in range(len(p)):
                                    if p[i][1] >= MIN_CLF_SCORE:
                                        #continue
                                        pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + "\t" + str(positionsPredicted[i]+1) + "\t" + str(round(p[i][1],5)) + "\t" + str(round(predictionArray[i][3],3)))
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
                        p = clf.predict_proba(predictionArray)
                        for i in range(len(p)):
                            if p[i][1] >= MIN_CLF_SCORE:
                                #continue
                                pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + "\t" + str(positionsPredicted[i]+1) + "\t" + str(round(p[i][1],5)) + "\t" + str(round(predictionArray[i][3],3)))
                        positionsPredicted = [] #bp positions of rows in predictionArray
                        nPred = 0

                if nPred > 0: #After last partial-batch of reads, clf remaining feature vectors
                    p = clf.predict_proba(predictionArray[:nPred])
                    for i in range(nPred):
                        if p[i][1] >= MIN_CLF_SCORE:
                                #continue
                                pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + "\t" + str(positionsPredicted[i]+1) + "\t" + str(round(p[i][1],5)) + "\t" + str(round(predictionArray[i][3],3)))
                                