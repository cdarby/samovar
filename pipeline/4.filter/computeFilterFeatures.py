#cdarby@jhu.edu
#updated 02-01-18
from __future__ import print_function

import pysam,argparse,pymp,sys
import numpy as np

if sys.version_info > (3,0):
	import _pickle as cPickle
else:
	import cPickle


def siteFeatures(depth, bases, bqs, mapqs, strands, readpos, clip, ind, XH, XL, XV, XR, XB, ASXS, pair):

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
	MAF_h0 = 0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if XH[i] == 0 and bases[i] == minor])/nphased
	MAF_h1 = 0 if nphased == 0 else 1.0*sum([1 for i in range(depth) if XH[i] == 1 and bases[i] == minor])/nphased
	
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
	Cfracclip = 1.0*sum([1 if i>0 else 0 for i in np.compress(mask,clip)]) / nC
	Cfracind = 1.0*sum([1 if i>0 else 0 for i in np.compress(mask,ind)]) / nC
	#CfracXT = 1.0*np.compress(mask,XT).sum()/nC
	CavgASXS = 1.0*np.compress(mask,ASXS).sum()/nC
	Cavgpair = 1.0*np.compress(mask,pair).sum()/nC

	mask_inv = [False if (mask[i] or XH[i] == None or (bases[i] != major and bases[i] != minor)) else True for i in range(depth) ]
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
	Nfracclip = 1.0*sum([1 if i>0 else 0 for i in np.compress(mask,clip)]) / nN
	Nfracind = 1.0*sum([1 if i>0 else 0 for i in np.compress(mask,ind)]) / nN
	#NfracXT = 1.0*np.compress(mask_inv,XT).sum()/nN
	NavgASXS = 1.0*np.compress(mask_inv,ASXS).sum()/nN
	Navgpair = 1.0*np.compress(mask_inv,pair).sum()/nN

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
	Mfracclip = 1.0*sum([1 if i>0 else 0 for i in np.compress(mask_M,clip)]) / nminor
	Mfracind = 1.0*sum([1 if i>0 else 0 for i in np.compress(mask_M,ind)]) / nminor
	#MfracXT = 1.0*np.compress(mask_M,XT).sum()/nminor
	MavgASXS = 1.0*np.compress(mask_M,ASXS).sum()/nminor
	Mavgpair = 1.0*np.compress(mask_M,pair).sum()/nminor

	mask_J = [True if bases[i] == major else False for i in range(depth)]
	Javgbq = 1.0*np.compress(mask_J,bqs).sum()/nmajor
	Javgmapq = 1.0*np.compress(mask_J,mapqs).sum()/nmajor
	Jfracstrand = 1.0*np.compress(mask_J,strands).sum()/nmajor
	Javgpos = 1.0*np.compress(mask_J,readpos).sum()/nmajor
	Javgclip = 1.0*np.compress(mask_J,clip).sum()/nmajor
	Javgind = 1.0*np.compress(mask_J,ind).sum()/nmajor
	Jfracclip = 1.0*sum([1 if i>0 else 0 for i in np.compress(mask_J,clip)]) / nmajor
	Jfracind = 1.0*sum([1 if i>0 else 0 for i in np.compress(mask_J,ind)]) / nmajor
	#JfracXT = 1.0*np.compress(mask_J,XT).sum()/nmajor
	JavgASXS = 1.0*np.compress(mask_J,ASXS).sum()/nmajor
	Javgpair = 1.0*np.compress(mask_J,pair).sum()/nmajor

	weightedMbq = 1.0*np.compress(mask_M,bqs).sum()/(Mavgbq * nminor + Javgbq * nmajor)
	Mavgibq = 1-weightedMbq #Should rm this feature - redundant. 
	
	features = [depth,fracphased,fracstrand,frach,MAF,MAF_phased,
		nC,fracC,Cavgbq,Cfrach,Cavgmapq,Cfracstrand,CMAF,CavgXL,CavgXV,CavgXR,CavgXB,Cavgpos,Cavgclip,Cavgind,Cfracclip,Cfracind,CavgASXS,Cavgpair,
		Navgbq,Nfrach,Navgmapq,Nfracstrand,NMAF,NavgXL,NavgXV,NavgXR,NavgXB,Navgpos,Navgclip,Nfracclip,Nfracind,Navgind,NavgASXS,Navgpair,
		weightedCbq,Cavgibq,
		Mavgbq,Mavgmapq,Mfracstrand,Mavgpos,Mavgclip,Mavgind,Mfracclip,Mfracind,Mavgibq,weightedMbq,MavgASXS,Mavgpair,
		Javgbq,Javgmapq,Jfracstrand,Javgpos,Javgclip,Javgind,Jfracclip,Jfracind,JavgASXS,Javgpair,MAFnorm,fracOtherAllele]

	return features

def readIntervalFile(fname):
	R = []
	with open(fname) as bedfile:
		region = bedfile.readline().strip().split()
		while len(region) > 0:
			R.append(region)
			region = bedfile.readline().strip().split()
	return R


NT = "ACGTN"

parser = argparse.ArgumentParser(description='python siteFeatures.py --bam diploid_output/phased_possorted_bam.bam --bed haploidPos.bed --out diploidBam_haploidPos.pileup')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name (output from classifySite)',required=True)
parser.add_argument('--nproc',help="parallelism",required=False,default=3)
parser.add_argument('--featuredict', help='.pkl file from calcLinkedReadFeatures.py' ,required=False,default=None)

#parser.add_argument('--out', help='Output file name',required=True)
args = parser.parse_args()
nproc = int(args.nproc)

R = readIntervalFile(args.bed)

if args.featuredict is not None:
	featureDict = cPickle.load(open(args.featuredict,'rb'))
	NROWS = max(featureDict.keys()) + 1
	featureArr = pymp.shared.array((NROWS,4),dtype='float32') #default float64
	for (k,v) in featureDict.items():
		featureArr[k][0] = v[0]
		featureArr[k][1] = v[1]
		featureArr[k][2] = v[2]
		featureArr[k][3] = v[3]
	featureDict = None #so that it's not copied by each parallel
else:
	featureArr = None

#outfile = open(args.out + args.bed +".filteringfeatures", "w")

with pymp.Parallel(nproc,if_=(nproc > 1)) as pymp_context:
	samfile = pysam.AlignmentFile(args.bam, "rb")
	for region_index in pymp_context.xrange(len(R)): #pymp.xrange returns an iterator and corresponds to dynamic scheduling.
		#pymp_context.print(pymp_context.thread_num)
		region = R[region_index] 
		chrom = region[0]
		start = int(region[1])# -1 # -1 for MosaicHunter
		end = int(region[2])
		clfscore = region[3]
		#print region
		for pileupcolumn in samfile.pileup(chrom, start, end, stepper="all"): 
			#allreads which overlap the region are returned. The first base returned will be the first base of the first read not necessarily the first base of the region used in the query.
			#skip reads in which any of the following flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
			if pileupcolumn.pos < start:
				continue
			if pileupcolumn.pos >= end:
				break

			readpos = []
			bases = []
			bqs = []
			mapqs = []
			clip = []
			ind = []
			strands = []

			XH = []
			XL = []
			XV = []
			XR = []
			XB = []
			#XT = []
			ASXS = []
			pair = []

			depth = 0
			for pileupread in pileupcolumn.pileups:
				aln = pileupread.alignment

				#base and base quality
				if pileupread.is_del and not pileupread.is_refskip:
					continue     
				else:
					depth += 1
					bases.append(str(aln.query_sequence[pileupread.query_position]))
				bqs.append(aln.query_qualities[pileupread.query_position])
				mapqs.append(aln.mapping_quality)
				readpos.append(min((aln.query_alignment_end - pileupread.query_position_or_next),(pileupread.query_position_or_next - aln.query_alignment_start)))
				clip.append(aln.query_length - len(aln.query_alignment_sequence))

				if aln.is_reverse: strands.append(0)
				else: strands.append(1)

				if aln.is_proper_pair: pair.append(0)
				else: pair.append(1)
				
				n_ind = 0
				for (op,ln) in aln.cigartuples:
					if op == 1 or op == 2 or op == 3:
						n_ind += ln
				ind.append(n_ind)
				
				#XT.append(aln.get_tag("XT"))
				ASXS.append(aln.get_tag("AS") - aln.get_tag("XS"))
				
				# linked read-quality features
				if aln.has_tag("HP"): #Longranger read phasing tag
					XH.append(aln.get_tag("HP")-1)
				else: 
					XH.append(None)
				if featureArr is not None and aln.has_tag("MI"): 
					MI = aln.get_tag("MI")
					try:
						M = featureArr[MI] 
						if M[0] < 1: #initialized to 0
							XR.append(0)
							XB.append(0)
							XL.append(0) #or None? Change in calcLinkedReadFeatures
							XV.append(0)
						else:
							XR.append(M[0]) #nreads
							XB.append(M[1]) #span
							XL.append(M[2]) #read agreement fraction
							XV.append(M[3]) #nvariants
					except IndexError:
						XR.append(0)
						XB.append(0)
						XL.append(0) #or None? Change in calcLinkedReadFeatures
						XV.append(0)
				else: 
					XR.append(0)
					XB.append(0)
					XL.append(0) #or None? Change in calcLinkedReadFeatures
					XV.append(0)

			if depth <= 0: 
				continue

			featurevector = siteFeatures(depth, bases, bqs, mapqs, strands, readpos, clip, ind, XH, XL, XV, XR, XB, ASXS, pair)
			if featurevector is not None:
				pymp_context.print(chrom + "\t" + str(start) + "\t" + clfscore + "\t" + "\t".join([str(x) for x in featurevector]))
