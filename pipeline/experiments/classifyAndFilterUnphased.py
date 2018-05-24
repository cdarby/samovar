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
#MIN_FRAC_PHASED = 0.5
MIN_FRAC_HAP1 = 0.3
MAX_FRAC_OTHERALLELE = 0.05
MIN_MAF = 0.05
MIN_HAPDISCORD_READS = 4
MIN_HAPDISCORD_READS_HAP = 0.1
MIN_AVG_POS_FROM_END = 10 #not for haplotype-concordant reads
MAX_AVG_INDELS = 5
MIN_CLF_SCORE = 0.75

SOFT_CLIP = True
GOOD_PAIR = True
STRAND = True

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
    if read.is_proper_pair: pair = 0
    else: pair = 1
    if read.is_reverse: strands = 0
    else: strands = 1      
    
    return (clip,clipLen,ind,pair,strands) 

def siteFeatures(S):
    (bases, bqs, readpos, clip, ind, pair, strands, clipLen) = (S.bases, S.bqs, S.readpos, S.clip, S.ind, S.pair, S.strands, S.clipLen)

    depth = len(bases)
    if depth < MIN_DEPTH: 
        #pymp_context.print("Depth")
        return None

    #features for whole site
    
    values, counts = np.unique(bases,return_counts=True)
    if len(values) == 1: return None #no minor allele bases
    basecounts = zip(values,counts)
    bases_ranked = sorted(basecounts, key=lambda x:x[1])
    major = bases_ranked[-1][0] #most frequent nt
    minor = bases_ranked[-2][0] #if bases_ranked[-2][0] != "N" else bases_ranked[-3][0]
    nmajor = bases_ranked[-1][1]
    nminor = bases_ranked[-2][1] #if bases_ranked[-2][0] != "N" else bases_ranked[-3][1]
    MAF = 1.0*nminor/depth
    if MAF < MIN_MAF: return None
    fracOtherAllele = 1.0*(depth - nmajor - nminor)/depth
    if fracOtherAllele > MAX_FRAC_OTHERALLELE: return None

    #non-linked features - "M-reads" are those with minor allele (like C-reads), "J-reads" are those with major allele (like N-reads)
    mask_M = [True if bases[i] == minor else False for i in range(depth)]
    Mavgpos = 1.0*np.compress(mask_M,readpos).sum()/nminor

    Mclip = np.compress(mask_M,clip).sum()
    Mpair = np.compress(mask_M,pair).sum()
    Mstrand = np.compress(mask_M,strands).sum()
    if Mavgpos < MIN_AVG_POS_FROM_END or (SOFT_CLIP and Mclip == nminor) or (GOOD_PAIR and Mpair == nminor) or (STRAND and (Mstrand == 0 or Mstrand == nminor)): #or Mavgind > MAX_AVG_INDELS 
        #pymp_context.print("M-read filter")
        return None

    mask_J = [True if bases[i] == major else False for i in range(depth)]
    Javgpos = 1.0*np.compress(mask_J,readpos).sum()/nmajor

    Jclip = np.compress(mask_J,clip).sum()
    Jpair = np.compress(mask_J,pair).sum()
    Jstrand = np.compress(mask_J,strands).sum()
    if Javgpos < MIN_AVG_POS_FROM_END or (SOFT_CLIP and Jclip == nmajor) or (GOOD_PAIR and Jpair == nmajor) or (STRAND and (Jstrand == 0 or Jstrand == nmajor)): #or Javgind > MAX_AVG_INDELS 
        #pymp_context.print("J-read filter")
        return None

    #if passes filters, do these computations
    Mavgbq = 1.0*np.compress(mask_M,bqs).sum()/nminor
    Mavgclip = 1.0*np.compress(mask_M,clipLen).sum()/nminor
    Mavgind = 1.0*np.compress(mask_M,ind).sum()/nminor

    Javgbq = 1.0*np.compress(mask_J,bqs).sum()/nmajor
    Javgclip = 1.0*np.compress(mask_J,clipLen).sum()/nmajor
    weightedMbq = 1.0*np.compress(mask_M,bqs).sum()/(Mavgbq * nminor + Javgbq * nmajor)
    Javgind = 1.0*np.compress(mask_J,ind).sum()/nmajor

    #features for classifier
    features = np.array([depth,MAF,Mavgbq,Mavgpos,Mavgclip,Mavgind,
        weightedMbq,Javgbq,Javgpos,Javgclip,Javgind])
    return features

parser = argparse.ArgumentParser(description='python classify.py --bam sample.bam --bed genome.bed --clf clf.pkl --nproc 8')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name',required=True)
parser.add_argument('--clf',help="random forest scikit-learn classifier",required=True)
parser.add_argument('--nproc',help="parallelism",required=False,default=3)
#parser.add_argument('--printproba',help="Print lines with at least this decision tree probability",required=False,default=0.9)

args = parser.parse_args()
Site = namedtuple('Site','pos bases bqs readpos clip clipLen ind HP pair strands')


NPAR = int(args.nproc)
READBATCH = 100000
CLFBATCH = 1000
INTERVALBATCH = 10000
MIN_CLF_SCORE = MIN_CLF_SCORE
NFEAT = 11

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
            for region in pymp_context.iterate(R): #pymp.xrange returns an iterator and corresponds to dynamic scheduling.
                #region = R[region_index]     
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
                    (clip,clipLen,ind,pair,strands) = readFeatures(read)

                    aligned_pairs = read.get_aligned_pairs()
                    for p in aligned_pairs:
                        if p[0] == None or p[1] == None:# or read.query_sequence[p[0]] == "N":
                            continue
                        querypos = p[0]
                        refpos = p[1]

                        if refpos in SITES:
                            S = SITES[refpos]
                        else:
                            S = Site(refpos,[],[],[],[],[],[],[],[],[])
                            SITES[refpos] = S

                        S.bases.append(read.query_sequence[querypos])
                        S.bqs.append(read.query_qualities[querypos])
                        S.readpos.append(min((read.query_alignment_end - querypos),(querypos - read.query_alignment_start)))
                        S.clip.append(clip)
                        S.ind.append(ind)
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
                                        pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + "\t" + str(positionsPredicted[i]+1) + "\t" + str(round(p[i][1],5)) + "\t" + str(round(predictionArray[i][1],3)))
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
                                pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + "\t" + str(positionsPredicted[i]+1) + "\t" + str(round(p[i][1],5)) + "\t" + str(round(predictionArray[i][1],3)))
                        positionsPredicted = [] #bp positions of rows in predictionArray
                        nPred = 0

                if nPred > 0: #After last partial-batch of reads, clf remaining feature vectors
                    p = clf.predict_proba(predictionArray[:nPred])
                    for i in range(nPred):
                        if p[i][1] >= MIN_CLF_SCORE:
                                #continue
                                pymp_context.print(chrom + "\t" + str(positionsPredicted[i]) + "\t" + str(positionsPredicted[i]+1) + "\t" + str(round(p[i][1],5)) + "\t" + str(round(predictionArray[i][1],3)))
                                