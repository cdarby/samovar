#cdarby@jhu.edu
#updated 02-01-18
from __future__ import print_function

import pysam,argparse,sys,string,random
import numpy as np
from sklearn.cluster import KMeans
from scipy.stats import fisher_exact

if sys.version_info > (3,0):
	import _pickle as cPickle
else:
	import cPickle

def readIntervalFile(fname):
	R = []
	with open(fname) as bedfile:
		region = bedfile.readline().strip().split()
		while len(region) > 0:
			R.append(region)
			region = bedfile.readline().strip().split()
	return R

def readVCF(fname):
	V = set()
	with open(fname) as vcf:
		fields = vcf.readline().strip().split()
		while len(fields) > 0:
			if "#" in fields[0]: 
				fields = vcf.readline().strip().split()
				continue
			V.add(fields[0] + ":" + str(int(fields[1])-1)) #Convert to 0-indexed
			fields = vcf.readline().strip().split()
	return V


NT = "ACGTN"

parser = argparse.ArgumentParser(description='python siteFeatures.py --bam diploid_output/phased_possorted_bam.bam --bed haploidPos.bed --out diploidBam_haploidPos.pileup')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name (output from classifySite)',required=True)
parser.add_argument('--ref',help="reference genome",required=True)
parser.add_argument('--vcfavoid',help="VCF of variants to NOT use as linkage",required=True)

#parser.add_argument('--featuredict', help='.pkl file from calcLinkedReadFeatures.py' ,required=True)

#parser.add_argument('--out', help='Output file name',required=True)
args = parser.parse_args()


R = readIntervalFile(args.bed)
V = readVCF(args.vcfavoid)
'''
featureDict = cPickle.load(open(args.featuredict,'rb'))
NROWS = max(featureDict.keys()) + 1
featureArr = pymp.shared.array((NROWS,4),dtype='float32') #default float64
for (k,v) in featureDict.items():
	featureArr[k][0] = v[0]
	featureArr[k][1] = v[1]
	featureArr[k][2] = v[2]
	featureArr[k][3] = v[3]
featureDict = None #so that it's not copied by each parallel
'''
reference = pysam.FastaFile(args.ref)
samfile = pysam.AlignmentFile(args.bam, "rb")
#for region_index in pymp_context.xrange(len(R)): #pymp.xrange returns an iterator and corresponds to dynamic scheduling.
for region in R:
	#region = R[region_index] 
	chrom = region[0]
	start = int(region[1])
	#end = int(region[2])
	alignments = []
	alignmentStarts = []
	alignmentEnds = []
	readNames = []
	basesAtSite = []
	pileupStart = None
	pileupEnd = None
	XH = [] 

	#print region
	for pileupcolumn in samfile.pileup(chrom, start, start+1, stepper="all"): 
		#if pileupStart == None: pileupStart = pileupcolumn.pos
		#allreads which overlap the region are returned. The first base returned will be the first base of the first read not necessarily the first base of the region used in the query.
		#Does this include soft clip on left side??
		#skip reads in which any of the following flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
		if pileupcolumn.pos < start:
			continue
		if pileupcolumn.pos > start:
			break

		for pileupread in pileupcolumn.pileups:
			aln = pileupread.alignment

			#base and base quality
			if pileupread.is_del and not pileupread.is_refskip:
				continue     #should do something with indels at the site in question?
			
			
			if aln.has_tag("HP"): #Longranger read phasing tag
				XH.append(aln.get_tag("HP")-1)
			else: 
				XH.append(None)

			readNames.append(aln.query_name)
			qstart = aln.query_alignment_start	
			rstart = aln.reference_start
			if pileupStart is None or rstart < pileupStart: pileupStart = rstart
			qend = aln.query_alignment_end
			rend = aln.reference_end		
			if pileupEnd is None or rend > pileupEnd: pileupEnd = rend
			qSeqAln = aln.query_alignment_sequence
			qSeqUnaln = aln.query_sequence
			#alignment location
			
			#thisReadFeatures.append(min((qend - pileupread.query_position_or_next),(pileupread.query_position_or_next - qstart))) #distance of base from end of aln
			#thisReadFeatures.append(qend) #aln end
			#thisReadFeatures.append(qstart) #aln start
			#thisReadFeatures.append(len(qSeqUnaln) - len(qSeqAln)) # num soft-clip
			#readsAtSite.append(thisReadFeatures)
			
			#mismatches and indels
			#How to encode for feature vector? Need entries from first pileupcolumn.pos to last pileupcolumn.pos
			alignmentStarts.append(rstart)
			alignmentEnds.append(rend)
			refSeq = reference.fetch(reference=chrom,start=rstart,end=rend).upper()
			refIdx = 0
			thisReadAln = [] #corresponds to positions in reference
			P = aln.get_aligned_pairs()
			for (q,r) in P: #position in query, position in reference
				if (q is not None and q >= len(qSeqUnaln)) or refIdx >= len(refSeq): break #should refIdx get too long?
				if r == start: #don't include the site of interest - see if they cluster otherwise.
					if q == None or qSeqUnaln[q] != refSeq[refIdx]: basesAtSite.append(1)
					else: basesAtSite.append(0)
					refIdx += 1
					continue
				elif r == None:
					if len(thisReadAln) > 0: thisReadAln[-1] = 1
					#insertion in query - change previous position in ref.
				elif q == None:
					thisReadAln.append(1)
					#deletion in query
					refIdx += 1
				elif qSeqUnaln[q] != refSeq[refIdx]:
					thisReadAln.append(1)
					#mismatch
					refIdx += 1
				else:
					thisReadAln.append(0)
					refIdx += 1
			alignments.append(thisReadAln)


	#identify c-reads
	depth = len(XH)
	if depth <= 0: continue
	if len(XH) != len(basesAtSite):
		print("\t".join([chrom, "None", "1", "1", "None"] + region))
		continue 
	change_to_hom_mask = [True if (XH[i] != None and basesAtSite[i] == 0) else False for i in range(depth)] #ignore non major/minor bases
	change_to_het_left_mask = [True if (XH[i] == 1 and basesAtSite[i] == 1) or (XH[i] == 0 and basesAtSite[i] == 0) else False for i in range(depth) ]
	change_to_het_right_mask = [True if (XH[i] == 0 and basesAtSite[i] == 1) or (XH[i] == 1 and basesAtSite[i] == 0) else False for i in range(depth) ]
	change_to_hom = sum(change_to_hom_mask)
	change_to_het_right = sum(change_to_het_right_mask)
	change_to_het_left = sum(change_to_het_left_mask)
	if change_to_hom <= change_to_het_left and change_to_hom <= change_to_het_right:
		mask = change_to_hom_mask
	elif change_to_het_left <= change_to_hom  and change_to_het_left <= change_to_het_right:
		mask = change_to_het_left_mask
	else:   
		mask = change_to_het_right_mask

	#pad all vectors with 0s to align the columns wrt reference
	# "ref" and "alt" are c/non-c reads
	totalLen = pileupEnd - pileupStart + 2
	depthAlt = [0 for _ in range(totalLen)]
	depthRef = [0 for _ in range(totalLen)]
	sumRef = [0 for i in range(totalLen)] 
	sumAlt = [0 for i in range(totalLen)]
	readsAtSite = [None for _ in range(depth)]

	for (i,seq) in enumerate(alignments):
		L = len(seq)
		readsAtSite[i] = [0 for _ in range(alignmentStarts[i] - pileupStart)] + [1] + seq + [1] + [0 for _ in range(pileupEnd - alignmentStarts[i] - L )] #the 1's surrounding seq indicate the read ends - hopefully columns are still in line.
		if mask[i] == 0:
			for j in range(-1,L+1): 
				depthRef[j + alignmentStarts[i] - pileupStart + 1] += 1
			for j in range(totalLen):
				sumRef[j] += readsAtSite[i][j]
		else:
			for j in range(-1,L+1): 
				depthAlt[j + alignmentStarts[i] - pileupStart + 1] += 1
			for j in range(totalLen):
				sumAlt[j] += readsAtSite[i][j]

	# should cache fisher pvals but will calc each time for now
	min_pval = 1
	min_table = None
	min_start = 0
	for k in range(totalLen):
		if chrom + ":" + str(pileupStart+k) in V: continue #Called variant from VCF
		T = [[sumAlt[k],depthAlt[k]-sumAlt[k]],[sumRef[k],depthRef[k]-sumRef[k]]]
		(oddsratio, pvalue) = fisher_exact(T)
		if pvalue < min_pval:
			min_table = T
			min_pval = pvalue
			min_start = pileupStart+k
	print("\t".join([str(x) for x in [chrom, min_start, min_pval, min_pval*totalLen, min_table] + region ]) )


