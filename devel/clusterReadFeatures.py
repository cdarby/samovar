#cdarby@jhu.edu
#updated 02-01-18
from __future__ import print_function

import pysam,argparse,sys,string,random
import numpy as np
from sklearn.cluster import KMeans

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


NT = "ACGTN"

parser = argparse.ArgumentParser(description='python siteFeatures.py --bam diploid_output/phased_possorted_bam.bam --bed haploidPos.bed --out diploidBam_haploidPos.pileup')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name (output from classifySite)',required=True)
parser.add_argument('--ref',help="reference genome",required=True)
#parser.add_argument('--featuredict', help='.pkl file from calcLinkedReadFeatures.py' ,required=True)

#parser.add_argument('--out', help='Output file name',required=True)
args = parser.parse_args()


R = readIntervalFile(args.bed)
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
	end = int(region[2])
	readsAtSite = []
	alignments = []
	alignmentStarts = []
	readNames = []
	basesAtSite = []
	pileupStart = None
	pileupEnd = None 

	#print region
	for pileupcolumn in samfile.pileup(chrom, start, end, stepper="all"): 
		#if pileupStart == None: pileupStart = pileupcolumn.pos
		#allreads which overlap the region are returned. The first base returned will be the first base of the first read not necessarily the first base of the region used in the query.
		#Does this include soft clip on left side??
		#skip reads in which any of the following flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
		if pileupcolumn.pos < start:
			continue
		if pileupcolumn.pos >= end:
			break

		for pileupread in pileupcolumn.pileups:
			aln = pileupread.alignment

			#base and base quality
			if pileupread.is_del and not pileupread.is_refskip:
				continue     #should do something with indels at the site in question?
			
			thisReadFeatures = []	

			#mapq related
			thisReadFeatures.append(aln.mapping_quality)
			thisReadFeatures.append(aln.get_tag("XT"))
			thisReadFeatures.append(aln.get_tag("AS") - aln.get_tag("XS"))

			#linked read related
			if aln.has_tag("HP"): #Longranger read phasing tag
				thisReadFeatures.append(aln.get_tag("HP")-1)
			else: 
				continue
			'''	
			if aln.has_tag("MI"): 
				MI = aln.get_tag("MI")
				try:
					M = featureArr[MI] 
					if M[0] < 1: #initialized to 0
						thisReadFeatures.append(0)
						thisReadFeatures.append(0)
						thisReadFeatures.append(0) #or None? Change in calcLinkedReadFeatures
						thisReadFeatures.append(0)
					else:
						thisReadFeatures.append(M[0]) #nreads
						thisReadFeatures.append(M[1]) #span
						thisReadFeatures.append(M[2]) #read agreement fraction
						thisReadFeatures.append(M[3]) #nvariants
				except IndexError:
					thisReadFeatures.append(0)
					thisReadFeatures.append(0)
					thisReadFeatures.append(0) #or None? Change in calcLinkedReadFeatures
					thisReadFeatures.append(0)
			else: 
				thisReadFeatures.append(0)
				thisReadFeatures.append(0)
				thisReadFeatures.append(0) #or None? Change in calcLinkedReadFeatures
				thisReadFeatures.append(0)
			'''
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
			thisReadFeatures.append(min((qend - pileupread.query_position_or_next),(pileupread.query_position_or_next - qstart)))
			thisReadFeatures.append(qend)
			thisReadFeatures.append(qstart)
			thisReadFeatures.append(len(qSeqUnaln) - len(qSeqAln))
			readsAtSite.append(thisReadFeatures)
			
			#mismatches and indels
			#How to encode for feature vector? Need entries from first pileupcolumn.pos to last pileupcolumn.pos
			alignmentStarts.append(rstart)
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

	#pad all vectors with 0s to align the columns wrt reference
	totalLen = pileupEnd - pileupStart + 2
	depthAlt = [0 for _ in range(totalLen)]
	depthRef = [0 for _ in range(totalLen)]
	sumRef = [0 for i in range(totalLen)]
	sumAlt = [0 for i in range(totalLen)]

	#otherFeats = [x[:] for x in readsAtSite]
	for (i,seq) in enumerate(alignments):
		L = len(seq)
		#readsAtSite[i] = [0 for _ in range(alignmentStarts[i] - pileupStart)] + seq + [0 for _ in range(pileupEnd - alignmentStarts[i] - L )] + readsAtSite[i]
		readsAtSite[i] = [0 for _ in range(alignmentStarts[i] - pileupStart)] + [1] + seq + [1] + [0 for _ in range(pileupEnd - alignmentStarts[i] - L )] #the 1's surrounding seq indicate the read ends - hopefully columns are still in line.
		if basesAtSite[i] == 0:
			for j in range(-1,L+1): 
				depthRef[j + alignmentStarts[i] - pileupStart + 1] += 1
			for j in range(totalLen):
				sumRef[j] += readsAtSite[i][j]
		else:
			for j in range(-1,L+1): 
				depthAlt[j + alignmentStarts[i] - pileupStart + 1] += 1
			for j in range(totalLen):
				sumAlt[j] += readsAtSite[i][j]
		'''
		if basesAtSite[i] == 0:
			for j in range(L): 
				depthRef[j + alignmentStarts[i] - pileupStart + 1] += 1
			for j in range(totalLen):
				sumRef[j] += readsAtSite[i][j]
		else:
			for j in range(L): 
				depthAlt[j + alignmentStarts[i] - pileupStart + 1] += 1
			for j in range(totalLen):
				sumAlt[j] += readsAtSite[i][j]
		'''
	#for k in readsAtSite:
	#	print("".join([str(i) for i in k]))
	#print(pileupStart,pileupEnd,pileupEnd-pileupStart,totalLen,len(readsAtSite),sum(basesAtSite))
	
	'''
	consensusDiff = []
	#positions = []
	for i in range(totalLen):
		if sumRef[i] != sumAlt[i]:
			d = abs(sumRef[i] - sumAlt[i])
			consensusDiff.append(min(depthAlt[i]+depthRef[i]-d,d))
			#positions.append(i+pileupStart)
	#print(zip(positions,consensusDiff))
	print(consensusDiff)
	'''
	'''
	r = random.shuffle(readsAtSite)
	sumRand1 = [0 for i in range(totalLen)]
	sumRand2 = [0 for i in range(totalLen)]
	for i in range(sum(basesAtSite)):
		for j in range(len(sumRand1)): sumRand1[j] += readsAtSite[i][j]
	for i in range(sum(basesAtSite),len(readsAtSite)):
		for j in range(len(sumRand2)): sumRand2[j] += readsAtSite[i][j]
	
	randDiff = []
	#positions = []
	for i in range(totalLen):
		if sumRand1[i] != sumRand2[i]:
			d = abs(sumRand1[i] - sumRand2[i])
			randDiff.append(min(depthAlt[i]+depthRef[i]-d,d))
			#positions.append(i+pileupStart)
	#print(zip(positions,consensusDiff))
	print(randDiff)
	'''

	
	#clusterRef = [readsAtSite[i] + otherFeats[i] for i in range(len(readsAtSite)) if basesAtSite[i] == 0]
	#clusterAlt = [readsAtSite[i] + otherFeats[i] for i in range(len(readsAtSite)) if basesAtSite[i] == 1]
	#totalLen = len(clusterRef[0])
	avgRef = [1.0*sumRef[i]/depthRef[i] if depthRef[i] > 0 else 0 for i in range(totalLen)]
	avgAlt = [1.0*sumAlt[i]/depthAlt[i] if depthAlt[i] > 0 else 0 for i in range(totalLen)]
	avgDiff = [str(abs(avgAlt[i] - avgRef[i])) for i in range(totalLen)]
	#print(chrom + "\t" + str(start) + "\t" + "\t".join(avgDiff))
	#print(chrom + "\t" + str(start) + "\t" + "\t".join([str(x) for x in sumAlt]) + "\t" + "\t".join([str(x) for x in sumRef]))
	print(chrom + "\t" + str(start) + "\t" + "\t".join([str(1.0*sumAlt[i]/(sumAlt[i]+sumRef[i])) if (sumRef[i]+sumAlt[i]) > 0 else "0" for i in range(totalLen)]) + "\t" + "\t".join([str(1.0*sumRef[i]/(sumAlt[i]+sumRef[i])) if (sumRef[i]+sumAlt[i]) > 0 else "0" for i in range(totalLen)]))
	#,1.0*sumRef[i]/(sumRef[i]+sumAlt[i]))
	#print(zip(sumRef,sumAlt))
	
	'''
	r = random.shuffle(readsAtSite)
	depthAlt = [0 for _ in range(totalLen)]
	depthRef = [0 for _ in range(totalLen)]
	for (i,seq) in enumerate(readsAtSite[:sum(basesAtSite)]):
		aln_start = readsAtSite[i].index(1)
		for j in range(L): 
			depthAlt[j + aln_start + 1] += 1
			sumAlt[j] += readsAtSite[i][j]

	for (i,seq) in enumerate(alignments[sum(basesAtSite):]):
		L = len(seq)
		readsAtSite[i] = [0 for _ in range(alignmentStarts[i] - pileupStart)] + [1] + seq + [1] + [0 for _ in range(pileupEnd - alignmentStarts[i] - L )]
		aln_start = readsAtSite[i].index(1)
		for j in range(L): 
			depthRef[j + aln_start + 1] += 1
			sumRef[j] += readsAtSite[i][j]

	avgRef = [1.0*sumRef[i]/depthRef[i] if depthRef[i] > 0 else 0 for i in range(totalLen)]
	avgAlt = [1.0*sumAlt[i]/depthAlt[i] if depthAlt[i] > 0 else 0 for i in range(totalLen)]
	avgDiff = [abs(avgAlt[i] - avgRef[i]) for i in range(totalLen)]		
	print(avgDiff)
	'''

	'''
	k = KMeans(n_clusters=4).fit_predict(readsAtSite)
	for (i,n) in enumerate(k):
		if basesAtSite[i] == 0:
			print(n,readNames[i])
	print()
	for (i,n) in enumerate(k):
		if basesAtSite[i] == 1:
			print(n,readNames[i])
	print()
	#sys.exit()
	'''


