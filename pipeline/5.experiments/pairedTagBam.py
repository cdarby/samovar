#cdarby@jhu.edu
#updated 10-23-17

#Access only variant sites, phase fragments, then label each read with the "PH" (paired haplotype) tag which is either 1 or 2

import pysam
import argparse, sys

parser = argparse.ArgumentParser(description='python pairedTagBam.py --bam sample.bam --out paired.bam --vcf sample.vcf')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--out', help='Output bam file name',required=True)
parser.add_argument('--vcf', help='Input vcf file name (a .vcf file from the sample)',required=True)

args = parser.parse_args()

alnfile = pysam.AlignmentFile(args.bam, "rb")
outfile = pysam.AlignmentFile(args.out, "wb",template=alnfile)

if not alnfile.check_index: 
	alnfile.close()
	sys.exit("Index the input bam file first")

#Read vcf 
variants = dict()
with open(args.vcf) as vcffile:
	line = vcffile.readline()
	while (line != ""):
		fields = line.split()
		if "#" not in fields[0] and "PASS" in fields[6] and len(fields[3]) == 1 and len(fields[4]) == 1:
			if fields[9][:3] == "0|1":
				k = fields[0] + ":" + str(int(fields[1])-1)
				variants[k] = fields[3] + fields[4]
			elif fields[9][:3] == "1|0":
				k = fields[0] + ":" + str(int(fields[1])-1)
				variants[k] = fields[4] + fields[3]
		line = vcffile.readline()

print("Read VCF")
print(str(len(variants)))

namesH0 = set()
namesH1 = set()

nv = 0
for (v,bases) in variants.items():
	nv += 1
	if nv % 10000 == 0: print(str(nv))
	(chrom,pos) = v.split(":")
	pos = int(pos)
	for pileupcolumn in alnfile.pileup(chrom,pos,pos+1):
		#default stepper=all which does "skip reads in which any of the following flags are set: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP"
		if pileupcolumn.pos > pos:
			break
		elif pileupcolumn.pos != pos:
			continue
		for pileupread in pileupcolumn.pileups:
			if not (pileupread.is_del or pileupread.is_refskip):
				#cov += 1
				aln = pileupread.alignment
				base = aln.query_sequence[pileupread.query_position] #TODO verify this
				if base == bases[0]: #hap 0
					namesH0.add(aln.query_name)
				elif base == variants[v][1]: #hap 1
					namesH1.add(aln.query_name)

print("Made phasing dictionary")
print(str(len(namesH0)),str(len(namesH1)))

nwritten=0
nohap=0
conflict=0
#Label reads
for r in alnfile.fetch():
	if nwritten % 100000 == 0: print(str(nwritten))
	N = r.query_name
	if N in namesH0:
		if N in namesH1: #mates have conflicting hap.
			outfile.write(r)
			nwritten+=1
			conflict += 1
			continue
		else: #H0
			r.set_tag("PH",1)
			outfile.write(r)
			nwritten+=1
			continue
	elif N in namesH1: #H1
		r.set_tag("PH",2)
		outfile.write(r)
		nwritten+=1
		continue
	else: #Neither
		outfile.write(r)
		nwritten+=1
		nohap += 1
		continue

print(str(nwritten),str(conflict),str(nohap))

alnfile.close()
outfile.close()
