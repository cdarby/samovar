#cdarby@jhu.edu
#updated 01-24-18

#Python 2.7.13 (tested 01-24-18)
#Python 3.6.2 (tested 01-24-18)

from __future__ import print_function

import pysam
import argparse, sys

if sys.version_info > (3,0):
	import _pickle as cPickle
else:
	import cPickle

def readVariants(vcf, isLR):
	variants = dict()
		#chrom : [position,letterAt0,letterAt1]
	with open(vcf) as vcffile:
		line = vcffile.readline()
		while (line != ""):
			fields = line.split()
			if isLR == None:
				if len(fields) > 1 and fields[0] != "BLOCK:" and fields[8] != "1" and len(fields[5]) == 1 and len(fields[6]) == 1: #some are separator lines; also non-pruned; also a SNP
					if fields[3] in variants: 
						V = variants[fields[3]]
					else: 
						variants[fields[3]] = [] 
						V = variants[fields[3]]
					if fields[1] + fields[2] == "01":
						V.append((int(fields[4])-1,fields[5],fields[6]))
					else: 
						V.append((int(fields[4])-1,fields[6],fields[5]))
			else:
				if "#" not in fields[0] and len(fields[3]) == 1 and len(fields[4]) == 1 and "PASS" in fields[6]:
					if fields[0] in variants: 
						V = variants[fields[0]]
					else: 
						variants[fields[0]] = [] 
						V = variants[fields[0]]
					if fields[9][:3] == "0|1":
						V.append((int(fields[1])-1,fields[3],fields[4]))
					elif fields[9][:3] == "1|0":
						V.append((int(fields[1])-1,fields[4],fields[3]))
			line = vcffile.readline()
	return variants

#########################

parser = argparse.ArgumentParser(description='makeLinkedReadFeatureDict.py')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--out', help='Output .pkl data structure name',required=False,default="features.pkl")
parser.add_argument('--vcf', help='Input vcf file name (a .hap file if from HC2, a .vcf file if from longranger)',required=True)
parser.add_argument('--isLR', help='Put 1 if using the LR vcf, otherwise expects HC2 .hap file',required=False,default=None)

args = parser.parse_args()

alnfile = pysam.AlignmentFile(args.bam, "rb")

if not alnfile.check_index: 
	alnfile.close()
	sys.exit("Index the input bam file first")

variants = readVariants(args.vcf,args.isLR)
print("Read " + str(len(variants)) + " contigs from VCF")
chrom = None
V = None
moleculeObjects = dict()
molecules = dict()

leftmostVariantCovered = 0
nVars = 0
for aln in alnfile.fetch():
	#aln = read.alignment
	#if not aln.has_tag("HP"): continue #molecule not phased
	if aln.has_tag("MI"): MI = aln.get_tag("MI")
	else: continue
	if aln.is_duplicate or aln.is_qcfail: continue #or aln.is_secondary or aln.is_unmapped or aln.is_supplementary
	#N += 1
	#if N % 1000000 == 0: print(N)
	alnstart = aln.reference_start #Use multiple times - access once
	alnend = aln.reference_end #Use multiple times - access once

	if aln.reference_name != chrom:
		chrom = aln.reference_name
		if chrom not in variants: 
			print("Contig name " + chrom + " not found in VCF; stopping here")
			break
		V = variants[chrom]
		#print(alnstart)
		leftmostVariantCovered = 0
		nVars = len(V)
		while leftmostVariantCovered < nVars and V[leftmostVariantCovered][0] < alnstart:
			leftmostVariantCovered += 1
		#dump moleculeObjects into molecules
		for (MI,M) in moleculeObjects.items():
			if (M[3] + M[4]) > 0:
				agreement = 1.0*abs(M[3] - M[4])/(M[3] + M[4])  
			else : agreement = 0.0 #or None?
			molecules[MI] = (M[2], (M[1]-M[0]), agreement, len(M[5])) #nreads, dist, agreement
		moleculeObjects = dict()
		print("Starting contig " + chrom)
		print("Molecules so far " + str(len(molecules)))

	if MI in moleculeObjects:
		M = moleculeObjects[MI]
		M[2] += 1
		if alnend > M[1]: M[1] = alnend
	else:
		M = [alnstart,alnend,1,0,0,set()] #startpos, endpos, nreads, readsh0, readsh1
		moleculeObjects[MI] = M

	if leftmostVariantCovered >= nVars: 
		continue
	if V[leftmostVariantCovered][0] < alnstart: leftmostVariantCovered += 1
	v = leftmostVariantCovered
	if v >= nVars or V[v][0] >= alnend: continue
	currentPair = 0
	P = aln.get_aligned_pairs()
	nPairs = len(P)
	while currentPair < nPairs:
		if P[currentPair][0] == None or P[currentPair][1] == None:
			currentPair += 1
			continue
		if P[currentPair][1] == V[v][0]:
			base = aln.query_sequence[P[currentPair][0]]
			M[5].add(P[currentPair][1])
			if base == V[v][1]: #hap 0
				M[3] += 1
			elif base == V[v][2]: #hap 1
				M[4] += 1
			v += 1
			if v >= nVars or V[v][0] >= alnend: break
		currentPair += 1

alnfile.close()

for (MI,M) in moleculeObjects.items():
	if (M[3] + M[4]) > 0:
		agreement = 1.0*abs(M[3] - M[4])/(M[3] + M[4])  
	else : agreement = 0.0 #or None?
	molecules[MI] = (M[2], (M[1]-M[0]), agreement, len(M[5])) #nreads, dist, agreement

outfile = open(args.out, "wb")
cPickle.dump(molecules,outfile)
print("Data structure is " + str(sys.getsizeof(molecules)) + " bytes in memory")
print("Wrote " + str(len(molecules)) + " molecules")
outfile.close()
