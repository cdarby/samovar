#cdarby@jhu.edu
#updated 09-27-18

#Python 2.7.10 (tested 05-24-18)
#Python 3.4.2 (tested 05-24-18)
#pypy 6.0.0 linux_x86_64-portable / Python 2.7.13 (tested 05-24-18)

from __future__ import print_function

import simplesam,argparse,sys
import numpy as np

##### Filter Parameters #####
MIN_DEPTH = 16
MIN_FRAC_PHASED = 0.5
MIN_FRAC_HAP1 = 0.3
MAX_FRAC_OTHERALLELE = 0.05
MIN_MAF = 0.05
MIN_HAPDISCORD_READS = 4
MAX_HAPDISCORD_READS_HAP = 0.1
MIN_AVG_POS_FROM_END = 10 #filter is not used for haplotype-concordant reads

SOFT_CLIP = True
GOOD_PAIR = True
STRAND = True

##### Performance parameters ####
READBATCH = 50000
SITEBATCH = 50000
WINDOW_SIZE = 2000000

# Contigs from the reference genome that WILL be scanned 
CONTIGS = set(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"])

def readFeatures(read):
	ind = 0
	clipLen = 0
	for c in read.cigars:
		if c[1] == "S": clipLen += c[0]
		elif c[1] != "M": ind += c[0]
	if clipLen > 0: clip = 1
	else: clip = 0  
	if read.paired: pair = 0
	else: pair = 1
	if read.reverse: strands = 0
	else: strands = 1      
	
	return (clip,clipLen,ind,pair,strands)


def siteFeatures(S):
	[refpos, bases, bqs, readpos, clip, ind, pair, strands, clipLen] = S
	depth = len(bases)
	if depth < MIN_DEPTH: 
		return None

	#features for whole site
	
	values, counts = np.unique(bases,return_counts=True)
	if len(values) == 1: return None #no minor allele bases
	basecounts = zip(values,counts)
	bases_ranked = sorted(basecounts, key=lambda x:x[1])
	major = bases_ranked[-1][0] #most frequent nt
	minor = bases_ranked[-2][0] #second most frequent nt
	nmajor = bases_ranked[-1][1]
	nminor = bases_ranked[-2][1]
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
	if Mavgpos < MIN_AVG_POS_FROM_END or (SOFT_CLIP and Mclip == nminor) or (GOOD_PAIR and Mpair == nminor) or (STRAND and (Mstrand == 0 or Mstrand == nminor)):
		return None

	mask_J = [True if bases[i] == major else False for i in range(depth)]
	Javgpos = 1.0*np.compress(mask_J,readpos).sum()/nmajor
	Jclip = np.compress(mask_J,clip).sum()
	Jpair = np.compress(mask_J,pair).sum()
	Jstrand = np.compress(mask_J,strands).sum()
	if Javgpos < MIN_AVG_POS_FROM_END or (SOFT_CLIP and Jclip == nmajor) or (GOOD_PAIR and Jpair == nmajor) or (STRAND and (Jstrand == 0 or Jstrand == nmajor)): 
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
	features = [depth,MAF,Mavgbq,Mavgpos,Mavgclip,Mavgind,
        weightedMbq,Javgbq,Javgpos,Javgclip,Javgind]
	return features


def processRegion(R): #samtools format region string
	chrom = R.split(":")[0]
	global CONTIGS
	if chrom not in CONTIGS: return []
	start = int(R.split(":")[1].split("-")[0])
	end = int(R.split(":")[1].split("-")[1])
	predictionArray = [] #list of feature vectors, 0th element in each vector is chrom and 1st is bp position
	SITES = dict()
	last_site_eval = -1
	nreads = 0
	global fname
	with open(fname, 'rb') as filenameopen:
		samfile = simplesam.Reader(filenameopen,regions=R)
	   
		while True:
			try: #get next read
				read = samfile.next()
				
				if read.duplicate or not read.passing or read.secondary: continue
				nreads += 1
				(clip,clipLen,ind,pair,strands) = readFeatures(read)

				gappedSeq = read.gapped('seq') #one time compute
				gappedQual = read.gapped('qual')
				#Features for the whole read
				for (p,refpos) in enumerate(read.coords): # "genomic coordinates for the gapped alignment"
					if gappedSeq[p] == "-":
						continue
					refpos = p + read.pos

					if refpos in SITES:
						S = SITES[refpos]
					else:
						S = [refpos,[],[],[],[],[],[],[],[]] #bases, bqs, readpos, clip, ind, HP, ASXS, pair, strands, clipLen
						SITES[refpos] = S
					
					S[1].append(gappedSeq[p])
					S[2].append(ord(gappedQual[p])-33)
					S[3].append(min((len(gappedSeq)-p),p)) #not exact due to indel
					S[4].append(clip)
					S[5].append(ind)
					S[6].append(pair)
					S[7].append(strands)
					S[8].append(clipLen)
					
				if nreads >= READBATCH or len(SITES) > SITEBATCH: #batch size for evaluating sites where all reads seen -> memory savings
					for refpos in list(SITES): # This will produce a list from the keys of the dictionary that will not change during iteration
						if refpos >= max(last_site_eval+1,start) and refpos < min(read.pos,end):
							S = SITES[refpos]
						else: continue
						feature_vector = siteFeatures(S)               
						if type(feature_vector) != type(None):
							predictionArray.append([chrom, refpos] + feature_vector)
						del SITES[S[0]]
					last_site_eval = read.pos - 1
					nreads = 0
				
			except StopIteration: #no more reads
				
				for refpos in list(SITES): # This will produce a list from the keys of the dictionary that will not change during iteration
					if refpos >= max(last_site_eval+1,start) and refpos < min(read.pos,end):
						S = SITES[refpos]
					else: continue
					feature_vector = siteFeatures(S)               
					if type(feature_vector) != type(None):
						predictionArray.append([chrom, refpos-1] + feature_vector)
					del SITES[S[0]]
				break #all sites classified
		samfile.close()
		samfile.p.wait() #Prevent Z-status samtools processes
		sys.stderr.write(R + "\n")
		return predictionArray

parser = argparse.ArgumentParser(description='python filter.py --bam sample.bam --nproc 8')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--nproc',help="parallelism",required=False,default=1)

args = parser.parse_args()
NFEAT = 33
NPAR = int(args.nproc)
fname = args.bam

with open(fname, 'rb') as F:
	S = simplesam.Reader(F)
	regions = S.tile_genome(WINDOW_SIZE)
	S.close()
	S.p.wait()

if NPAR > 1:
	import multiprocessing as mp
	pool = mp.Pool(processes=NPAR,maxtasksperchild=1) #new process is started for each region, max NPAR active at a time
	for result in pool.imap_unordered(processRegion,regions):
		p = [str(x) for x in result] #predictionArray
		if len(p) > 0: print("\n".join(p))
else:
	for r in regions:
		p = [str(x) for x in processRegion(r)] #predictionArray
		if len(p) > 0: print("\n".join(p))
