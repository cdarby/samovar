#cdarby@jhu.edu
#updated 05-23-18

#Python 2.7.10 (tested 05-24-18)
#Python 3.4.2 (tested 05-24-18)
#pypy 6.0.0 linux_x86_64-portable / Python 2.7.13 (tested 05-24-18)

from __future__ import print_function

import simplesam,argparse,sys
import numpy as np
from fisher import pvalue
from pyfaidx import Fasta

if sys.version_info > (3,0):
	import _pickle as cPickle
else:
	import cPickle

MIN_FISHER_PVAL = 0.05

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


def processRegion(R): #bed interval -> samtools format region string
	global fname
	chrom = R[0]
	start = int(R[1])
	r = chrom + ":" + str(start+1) + "-" + str(start+2) #1-indexed
	alignments = []
	alignmentStarts = []
	alignmentEnds = []
	readNames = []
	basesAtSite = []
	pileupStart = None
	pileupEnd = None

	with open(fname, 'rb') as filenameopen:
		samfile = simplesam.Reader(filenameopen,regions=r)
		global fastafname
		reference = Fasta(fastafname) #Fasta() is apparently not "thread safe" to use as global and creating per-process does not affect performance
		numPhased = 0
		while True:
			try: #get next read
				read = samfile.next()
				if read.duplicate or not read.passing or read.secondary: continue
				readNames.append(read.qname)
				try:
					HP = read["HP"]
					numPhased += 1
				except:
					HP = None
				C = read.cigars
				if C[0][1] == "S":
					qstart = C[0][0]	
				else:
					qstart = 0
				rstart = read.pos
				if pileupStart is None or rstart < pileupStart: pileupStart = rstart
				if C[-1][1] == "S":
					qend = C[-1][0]	
				else:
					qend = 0				
				rend = rstart + len(read)
				if pileupEnd is None or rend > pileupEnd: pileupEnd = rend

				qSeqAln = read.gapped('seq')

				#mismatches and indels
				alignmentStarts.append(rstart)
				alignmentEnds.append(rend)

				refSeq = reference[R[0]][rstart-1:rend-1].seq.upper()
				thisReadAln = []
				P = read.coords #positions in reference for qSeqAln
				addedBase = False			
				try:
					for i in range(len(P)): #position in reference
						refPos = P[i]
						qChar = qSeqAln[i]
						if refPos == start: #don't include the site of interest - see if they cluster otherwise.
							if i > 0 and refPos-1 != P[i-1] or qChar != refSeq[refPos-rstart]: basesAtSite.append(1)
							else: basesAtSite.append(0)
							addedBase = True
						elif i > 0 and refPos-1 != P[i-1]:
							if len(thisReadAln) > 0: thisReadAln[-1] = 1
							#insertion in query - change previous position in ref.
						elif qChar == "-":
							thisReadAln.append(1)
							#deletion in query
						elif qChar != refSeq[refPos-rstart]:
							thisReadAln.append(1)
							#mismatch
						else:
							thisReadAln.append(0)
				except IndexError:
					sys.stderr.write(str(len(P))+"\n")
					sys.stderr.write(str(len(qSeqAln))+"\n")
					sys.stderr.write(str(len(refSeq))+"\n")
					sys.stderr.write(str(max(P))+"\n")
					sys.stderr.write(str(rstart)+"\n")
					sys.stderr.write(str(P)+"\n")
					sys.stderr.write(str(qSeqAln)+"\n")
					sys.stderr.write(str(refSeq)+"\n")
					return None

				alignments.append(thisReadAln)
				if not addedBase: basesAtSite.append(1) #read didn't align to that position
			except StopIteration: #no more reads
				depth = len(readNames)
				if depth <= 0: 
					predictionString = None
					sys.stderr.write(str(R)+"\n")
					break
				fracPhased = 1.0*numPhased/depth
				totalLen = pileupEnd - pileupStart + 2
				depthAlt = [0 for _ in range(totalLen)]
				depthRef = [0 for _ in range(totalLen)]
				sumRef = [0 for i in range(totalLen)] 
				sumAlt = [0 for i in range(totalLen)]
				readsAtSite = [None for _ in range(depth)]

				for (i,seq) in enumerate(alignments):
					L = len(seq)
					readsAtSite[i] = [0 for _ in range(alignmentStarts[i] - pileupStart)] + [1] + seq + [1] + [0 for _ in range(pileupEnd - alignmentStarts[i] - L )] #the 1's surrounding seq indicate the read ends - hopefully columns are still in line.
					if basesAtSite[i] == 0:
						for j in range(-1,L+1): 
							depthRef[j + alignmentStarts[i] - pileupStart + 1] += 1
						for j in range(totalLen):
							sumRef[j] += readsAtSite[i][j] #events on h0
					else:
						for j in range(-1,L+1): 
							depthAlt[j + alignmentStarts[i] - pileupStart + 1] += 1
						for j in range(totalLen):
							sumAlt[j] += readsAtSite[i][j] #events on h1

				min_pval = 1
				min_table = None
				min_start = 0
				global V
				for k in range(totalLen):
					if chrom + ":" + str(pileupStart+k) in V: continue #Called variant from VCF
					p = pvalue(sumAlt[k],depthAlt[k]-sumAlt[k],sumRef[k],depthRef[k]-sumRef[k])
					if p.two_tail < min_pval:
						min_pval = p.two_tail
				if min_pval < MIN_FISHER_PVAL: predictionString = None
				else: predictionString = ("\t".join(R + [str(min_pval),str(fracPhased)]))
				break

	samfile.close()
	samfile.p.wait() #Prevent Z-status samtools process
	return predictionString


parser = argparse.ArgumentParser(description='linkageFilter.py --bam example.bam --bed predictions.txt --vcfavoid example.vcf --ref genome.fa')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name (output from classify)',required=True)
parser.add_argument('--ref',help="reference genome",required=True)
parser.add_argument('--vcfavoid',help="VCF of variants to NOT use as linkage",required=False,default=None)
parser.add_argument('--nproc',help="parallelism",required=False,default=1)
args = parser.parse_args()

regions = readIntervalFile(args.bed)
if args.vcfavoid is not None:
	V = readVCF(args.vcfavoid)
else:
	V = set()

fastafname = args.ref
fname = args.bam
NPAR = int(args.nproc)

if NPAR > 1:
	import multiprocessing as mp
	pool = mp.Pool(processes=NPAR)
	for result in pool.imap_unordered(processRegion,regions):
		if result is not None: print(result)
else:
	for r in regions:
		result = processRegion(r)
		if result is not None: print(result)

