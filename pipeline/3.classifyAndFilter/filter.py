#cdarby@jhu.edu
#updated 05-23-18

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
	ASXS = read["AS"] - read["XS"]
	try:
		HP = int(read["HP"])-1
	except:
		HP = None  
	if read.paired: pair = 0
	else: pair = 1
	if read.reverse: strands = 0
	else: strands = 1      
	
	return (clip,clipLen,ind,HP,ASXS,pair,strands)


def siteFeatures(S):
	[refpos, bases, bqs, readpos, clip, ind, HP, ASXS, pair, strands, clipLen] = S
	depth = len(bases)
	if depth < MIN_DEPTH: 
		return None

	#features for whole site
	
	nh0 = sum([k == 0 for k in HP])
	nh1 = sum([k == 1 for k in HP])
	nphased = nh0 + nh1
	if nphased < MIN_HAPDISCORD_READS: 
		return None
	fracphased = 1.0*nphased/depth
	if fracphased < MIN_FRAC_PHASED: 
		return None
	frach0 = 1.0*nh0/nphased
	frach1 = 1.0*nh1/nphased
	frach = max(frach0,frach1)
	if frach < MIN_FRAC_HAP1 or (1-frach) < MIN_FRAC_HAP1: 
		return None
	
	values, counts = np.unique(bases,return_counts=True)
	if len(values) == 1: return None #no minor allele bases
	basecounts = zip(values,counts)
	bases_ranked = sorted(basecounts, key=lambda x:x[1])
	major = bases_ranked[-1][0] #most frequent nt
	minor = bases_ranked[-2][0] #second most frequent nt
	nmajor = bases_ranked[-1][1]
	nminor = bases_ranked[-2][1]
	MAF = 1.0*nminor/depth
	MAF_phased = 0.0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if bases[i] != major and HP[i] != None])/nphased
	if MAF < MIN_MAF: return None
	fracOtherAllele = 1.0*(depth - nmajor - nminor)/depth
	if fracOtherAllele > MAX_FRAC_OTHERALLELE: return None

	#linked features
	#haplotype discordant reads: need to be changed in the allele x haplotype matrix to make alleles on each haplotype uniform
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
	if not (Cfrach < MAX_HAPDISCORD_READS_HAP or (1-Cfrach) < MAX_HAPDISCORD_READS_HAP): 
		return None

	Cavgpos = 1.0*np.compress(mask,readpos).sum()/nC #average distance from end of alignments
	Cclip = np.compress(mask,clip).sum() #number of reads clipped
	Cpair = np.compress(mask,pair).sum() #number of reads with bad pair
	Cstrand = np.compress(mask,strands).sum() #number of reads on + strand
	if Cavgpos < MIN_AVG_POS_FROM_END or (SOFT_CLIP and Cclip == nC) or (GOOD_PAIR and Cpair == nC) or (STRAND and (Cstrand == 0 or Cstrand == nC)):
		return None

	mask_inv = [False if (mask[i] or HP[i] == None or (bases[i] != major and bases[i] != minor)) else True for i in range(depth) ]
	nN = sum(mask_inv)

	Navgpos = 1.0*np.compress(mask_inv,readpos).sum()/nN
	Npair = np.compress(mask_inv,pair).sum()
	Nclip = np.compress(mask_inv,clip).sum() #number of reads clipped
	Nstrand = np.compress(mask_inv,strands).sum() #number of reads on + strand
	if (SOFT_CLIP and Nclip == nN) or (GOOD_PAIR and Npair == nN) or (STRAND and (Nstrand == 0 or Nstrand == nN)): 
		return None

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

	C_bases = np.compress(mask,bases)
	values, counts = np.unique(C_bases,return_counts=True)
	basecounts = zip(values,counts)
	bases_ranked = sorted(basecounts, key=lambda x:x[1])
	CMAF = 1.0*bases_ranked[-1][1]/depth
	CavgASXS = 1.0*np.compress(mask,ASXS).sum()/nC
	Cavgclip = 1.0*np.compress(mask,clipLen).sum()/nC #average number of bases clipped
	fracC = 1.0*nC/(nN+nC)
	Cavgbq = 1.0*np.compress(mask,bqs).sum()/nC
	Cavgind = 1.0*np.compress(mask,ind).sum()/nC #average number of indels in alignment

	N_bases = np.compress(mask_inv,bases)
	values, counts = np.unique(N_bases,return_counts=True)
	basecounts = zip(values,counts)
	bases_ranked = sorted(basecounts, key=lambda x:x[1])
	NMAF = 1.0*bases_ranked[-1][1]/depth
	NavgASXS = 1.0*np.compress(mask_inv,ASXS).sum()/nN
	Navgclip = 1.0*np.compress(mask_inv,clipLen).sum()/nN
	Navgind = 1.0*np.compress(mask_inv,ind).sum()/nN
	
	Navgbq = 1.0*np.compress(mask_inv,bqs).sum()/nN
	N_haps = np.compress(mask_inv,HP)
	Nfrach0 = 1.0*sum([1 for i in range(nN) if N_haps[i] == 0])/nN
	Nfrach1 = 1.0*sum([1 for i in range(nN) if N_haps[i] == 1])/nN
	Nfrach = max(Nfrach0,Nfrach1)
	weightedCbq = 1.0*nC*Cavgbq/(nC*Cavgbq+nN*Navgbq)

	Mavgbq = 1.0*np.compress(mask_M,bqs).sum()/nminor
	Mavgclip = 1.0*np.compress(mask_M,clipLen).sum()/nminor
	MavgASXS = 1.0*np.compress(mask_M,ASXS).sum()/nminor
	Mavgind = 1.0*np.compress(mask_M,ind).sum()/nminor  

	Javgbq = 1.0*np.compress(mask_J,bqs).sum()/nmajor
	Javgclip = 1.0*np.compress(mask_J,clipLen).sum()/nmajor
	JavgASXS = 1.0*np.compress(mask_J,ASXS).sum()/nmajor
	weightedMbq = 1.0*np.compress(mask_M,bqs).sum()/(Mavgbq * nminor + Javgbq * nmajor)
	Javgind = 1.0*np.compress(mask_J,ind).sum()/nmajor

	#features for classifier
	features = [minor,depth,fracphased,frach,MAF,MAF_phased,
		nC,fracC,Cavgbq,Cfrach,CMAF,
		Cavgpos,Cavgclip,Cavgind,CavgASXS,Navgbq,
		Nfrach,NMAF,Navgpos,Navgclip,Navgind,
		NavgASXS,weightedCbq,Mavgbq,Mavgpos,Mavgclip,
		Mavgind,MavgASXS,weightedMbq,Javgbq,Javgpos,
		Javgclip,Javgind,JavgASXS]
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
				(clip,clipLen,ind,HP,ASXS,pair,strands) = readFeatures(read)

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
						S = [refpos,[],[],[],[],[],[],[],[],[],[]] #bases, bqs, readpos, clip, ind, HP, ASXS, pair, strands, clipLen
						SITES[refpos] = S
					
					S[1].append(gappedSeq[p])
					S[2].append(ord(gappedQual[p])-33)
					S[3].append(min((len(gappedSeq)-p),p)) #not exact due to indel
					S[4].append(clip)
					S[5].append(ind)
					S[6].append(HP)
					S[7].append(ASXS)
					S[8].append(pair)
					S[9].append(strands)
					S[10].append(clipLen)
					
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
