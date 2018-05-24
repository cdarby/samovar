#cdarby@jhu.edu
#updated 05-23-18

#Python 2.7.10 (tested 05-24-18)
#Python 3.4.2 (tested 05-24-18)
#pypy 6.0.0 linux_x86_64-portable / Python 2.7.13 (tested 05-24-18)

from __future__ import print_function

import simplesam,argparse,sys,random
import numpy as np

NT = "ACGT"

def readIntervalFile(fname,simMode,hetMode,homMode):
	R = []
	with open(fname) as varfile:
		region = varfile.readline().strip().split()
		while len(region) > 0:
			if "#" in region[0]:
				region = varfile.readline().strip().split()
				continue
			chrom = region[0]
			
			if simMode: 
				vaf = float(region[2])
				site = int(region[1]) #varfile is 0-indexed
			elif hetMode: 
				if len(region[3]) != 1 or len(region[4]) != 1 or "PASS" not in region[6] or ("0|1" not in region[9] and "1|0" not in region[9]):
					region = varfile.readline().strip().split()
					continue
				vaf = None
				site = int(region[1]) - 1 #vcf is 1-indexed
			elif homMode: 
				if len(region[3]) != 1 or len(region[4]) != 1 or "PASS" not in region[6] or ("0|0" not in region[9] and "1|1" not in region[9]):
					region = varfile.readline().strip().split()
					continue
				vaf = None
				site = int(region[1]) - 1 #vcf is 1-indexed
			R.append([chrom,site,vaf])
			region = varfile.readline().strip().split()
	return R

def readFeatures(read):
	ind = 0
	clipLen = 0
	for c in read.cigars:
		if c[1] == "S": clipLen += c[0]
		elif c[1] != "M": ind += c[0]
	ASXS = read["AS"] - read["XS"]
	try:
		HP = int(read["HP"])-1
	except:
		HP = None      
	return (clipLen,ind,HP,ASXS)


def siteFeatures(S,vaf,simMode):
	[refpos, bases, bqs, readpos, clip, ind, HP, ASXS] = S
	depth = len(bases)
	if depth == 0: 
		return None

	#features for whole site
	
	nh0 = sum([k == 0 for k in HP])
	nh1 = sum([k == 1 for k in HP])
	nphased = nh0 + nh1
	if nphased == 0: return None
	fracphased = 1.0*nphased/depth
	frach0 = 1.0*nh0/nphased
	frach1 = 1.0*nh1/nphased
	frach = max(frach0,frach1)
	
	values, counts = np.unique(bases,return_counts=True)
	basecounts = zip(values,counts)
	bases_ranked = sorted(basecounts, key=lambda x:x[1])
	major = bases_ranked[-1][0] #most frequent nt

	if simMode:
		if len(values) == 1: 
			minor = major
			while minor == major:
				minor = NT[random.randint(0,3)] #both bounds inclusive
		else:
			minor = bases_ranked[-2][0]
	
		new_bases = []
		for (i,b) in enumerate(bases):
			if b != major or HP[i] == 0: 
				new_bases.append(b)
			elif HP[i] is None:
				if random.random() < vaf/2: #assume non-haplotyped reads are distributed equally on h1 and h0
					new_bases.append(minor)
				else:
					new_bases.append(b)    
			else: # HP[i] == 1:
				if random.random() < vaf: 
					new_bases.append(minor)
				else:
					new_bases.append(b)
		bases = new_bases

		values, counts = np.unique(bases,return_counts=True)
		#if len(values) == 1: return None #no minor allele bases
		basecounts = zip(values,counts)
		bases_ranked = sorted(basecounts, key=lambda x:x[1])
		major = bases_ranked[-1][0] #most frequent nt

	if len(values) == 1: return None #no minor allele bases
	minor = bases_ranked[-2][0]
	nmajor = bases_ranked[-1][1]
	nminor = bases_ranked[-2][1]
	MAF = 1.0*nminor/depth
	fracOtherAllele = 1.0*(depth - nmajor - nminor)/depth
	MAF_phased = 0.0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if bases[i] != major and HP[i] != None])/nphased #MAF of phased reads
	MAF_h0 = 0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if HP[i] == 0 and bases[i] == minor])/nphased
	MAF_h1 = 0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if HP[i] == 1 and bases[i] == minor])/nphased

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
	if nC == 0 or nC == nphased: #No C-reads or no N-reads
		return None 

	Cavgbq = 1.0*np.compress(mask,bqs).sum()/nC
	C_haps = np.compress(mask,HP)
	Cfrach0 = 1.0*sum([1 for i in range(nC) if C_haps[i] == 0])/nC
	Cfrach1 = 1.0*sum([1 for i in range(nC) if C_haps[i] == 1])/nC
	Cfrach = max(Cfrach0,Cfrach1)
	
	C_bases = np.compress(mask,bases)
	values, counts = np.unique(C_bases,return_counts=True)
	basecounts = zip(values,counts)
	bases_ranked = sorted(basecounts, key=lambda x:x[1])
	CMAF = 1.0*bases_ranked[-1][1]/depth
	
	Cavgpos = 1.0*np.compress(mask,readpos).sum()/nC
	Cavgclip = 1.0*np.compress(mask,clip).sum()/nC
	Cavgind = 1.0*np.compress(mask,ind).sum()/nC
	CavgASXS = 1.0*np.compress(mask,ASXS).sum()/nC

	mask_inv = [False if (mask[i] or HP[i] == None or (bases[i] != major and bases[i] != minor)) else True for i in range(depth) ]
	nN = sum(mask_inv)
	Navgbq = 1.0*np.compress(mask_inv,bqs).sum()/nN
	N_haps = np.compress(mask_inv,HP)
	Nfrach0 = 1.0*sum([1 for i in range(nN) if N_haps[i] == 0])/nN
	Nfrach1 = 1.0*sum([1 for i in range(nN) if N_haps[i] == 1])/nN
	Nfrach = max(Nfrach0,Nfrach1)
		
	N_bases = np.compress(mask_inv,bases)
	values, counts = np.unique(N_bases,return_counts=True)
	basecounts = zip(values,counts)
	bases_ranked = sorted(basecounts, key=lambda x:x[1])
	NMAF = 1.0*bases_ranked[-1][1]/depth

	Navgpos = 1.0*np.compress(mask_inv,readpos).sum()/nN
	Navgclip = 1.0*np.compress(mask_inv,clip).sum()/nN
	Navgind = 1.0*np.compress(mask_inv,ind).sum()/nN
	NavgASXS = 1.0*np.compress(mask_inv,ASXS).sum()/nN

	fracC = 1.0*nC/(nN+nC)
	weightedCbq = 1.0*nC*Cavgbq/(nC*Cavgbq+nN*Navgbq)
	
	#non-linked features - "M-reads" are those with minor allele (like C-reads), "J-reads" are those with major allele (like N-reads)
	mask_M = [True if bases[i] == minor else False for i in range(depth)]
	Mavgbq = 1.0*np.compress(mask_M,bqs).sum()/nminor
	Mavgpos = 1.0*np.compress(mask_M,readpos).sum()/nminor
	Mavgclip = 1.0*np.compress(mask_M,clip).sum()/nminor
	Mavgind = 1.0*np.compress(mask_M,ind).sum()/nminor
	MavgASXS = 1.0*np.compress(mask_M,ASXS).sum()/nminor

	mask_J = [True if bases[i] == major else False for i in range(depth)]
	Javgbq = 1.0*np.compress(mask_J,bqs).sum()/nmajor
	Javgpos = 1.0*np.compress(mask_J,readpos).sum()/nmajor
	Javgclip = 1.0*np.compress(mask_J,clip).sum()/nmajor
	Javgind = 1.0*np.compress(mask_J,ind).sum()/nmajor
	JavgASXS = 1.0*np.compress(mask_J,ASXS).sum()/nmajor

	weightedMbq = 1.0*np.compress(mask_M,bqs).sum()/(Mavgbq * nminor + Javgbq * nmajor)

	features = [depth,fracphased,frach,MAF,MAF_phased,
		nC,fracC,Cavgbq,Cfrach,CMAF,
		Cavgpos,Cavgclip,Cavgind,CavgASXS,Navgbq,
		Nfrach,NMAF,Navgpos,Navgclip,Navgind,
		NavgASXS,weightedCbq,Mavgbq,Mavgpos,Mavgclip,
		Mavgind,MavgASXS,weightedMbq,Javgbq,Javgpos,
		Javgclip,Javgind,JavgASXS]
	return features

def processRegion(R):
	global fname
	r = R[0] + ":" + str(R[1]) + "-" + str(R[1]+1)
	with open(fname, 'rb') as filenameopen:
		samfile = simplesam.Reader(filenameopen,regions=r)
	
		S = [R[1],[],[],[],[],[],[],[]]
		while True:
			try: #get next read
				read = samfile.next()
				if read.duplicate or not read.passing or read.secondary: continue
				#indel at site in question?
				(clip,ind,HP,ASXS) = readFeatures(read)
				gappedSeq = read.gapped('seq')
				S[1].append(gappedSeq[R[1]-read.pos])
				S[2].append(ord(read.gapped('qual')[R[1]-read.pos])-33)
				S[3].append(min((len(gappedSeq)-R[1]),R[1])) #not exact due to indel
				S[4].append(clip)
				S[5].append(ind)
				S[6].append(HP)
				S[7].append(ASXS)

			except StopIteration: #no more reads  
				samfile.close()
				samfile.p.wait() #Prevent Z-status samtools processes
				global simMode
				return siteFeatures(S,R[2],simMode)


parser = argparse.ArgumentParser(description='python simulate.py --bam sample.bam --varfile out.varfile --simulate --nproc 8')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--varfile', help='Input varfile name (.varfile if --simulate True, otherwise .vcf)',required=True)
parser.add_argument('--simulate', help='Specify iff simulation mode',action='store_true')
parser.add_argument('--het', help='Specify iff vcf is given and want het sites',action='store_true')
parser.add_argument('--hom', help='Specify iff vcf is given and want hom sites',action='store_true')
parser.add_argument('--max', help='Maximum sites to (successfully) simulate',required=False,default=None)
parser.add_argument('--nproc',help="parallelism",required=False,default=1)
args = parser.parse_args()

NPAR = int(args.nproc)
simMode = args.simulate
hetMode = args.het
homMode = args.hom
maxdata = int(args.max) if args.max is not None else None
regions = readIntervalFile(args.varfile,simMode,hetMode,homMode)
fname = args.bam

if maxdata is not None: #Limit on how many results to return as opposed to using the whole VCF file in --het and --hom when you only need a few thousand training examples
	i = 0
	if NPAR > 1:
		import multiprocessing as mp
		pool = mp.Pool(processes=NPAR)
		for result in pool.imap_unordered(processRegion,regions):
			if result is not None: 
				print("\t".join([str(x) for x in result]))
				i += 1
				if i > maxdata: sys.exit()
	else:
		for r in regions:
			result = processRegion(r)
			if result is not None: 
				print("\t".join([str(x) for x in result]))
				i += 1
				if i > maxdata: sys.exit()
else:
	if NPAR > 1:
		import multiprocessing as mp
		pool = mp.Pool(processes=NPAR)
		for result in pool.imap_unordered(processRegion,regions):
			if result is not None: print("\t".join([str(x) for x in result]))
	else:
		for r in regions:
			result = processRegion(r)
			if result is not None: print("\t".join([str(x) for x in result]))



