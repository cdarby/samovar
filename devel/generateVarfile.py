#cdarby@jhu.edu
#updated 01-24-18

#Python 2.7.13 (tested 01-24-18)
#Python 3.6.2 (tested 01-24-18)

from __future__ import print_function

import argparse,os,sys,random
from collections import namedtuple


parser = argparse.ArgumentParser(description='python generateVarfile.py --out out.varfile --vcf sample.vcf --fai genome.fa.fai --variantspacing 100000 --vafspacing 0.05')
parser.add_argument('--out', help='Output varfile name',required=False, default="out.varfile")
parser.add_argument('--vcf', help='Which sites to avoid (takes comma separated multiple file names)',required=False)
parser.add_argument('--fai', help='fasta file index of the genome to simulate',required=True)
parser.add_argument('--variantspacing', help='spacing between simulated sites',required=False,default=100000)
parser.add_argument('--vafspacing', help='spacing between variant allele fraction levels',required=False,default=0.05)

args = parser.parse_args()

avoid = set()
if args.vcf != None:
	A = args.vcf.split(",")
	for fname in A:
		with open(fname) as F:
			line = F.readline().strip().split()
			while len(line) != 0:
				if "#" not in line[0]:
					avoid.add(line[0] + ":" + str(int(line[1])-1)) #offset 1 for "avoid" sites because it's a VCF
				line = F.readline().strip().split()

vafspacing = float(args.vafspacing)
variantspacing = int(args.variantspacing)
with open(args.out,"w") as outfile, open(args.fai) as faifile:
	interval = faifile.readline().strip().split()
	while len(interval) > 0:
		seqname = interval[0]
		total_length = int(interval[1])
		pos = variantspacing
		vaf = vafspacing
		while pos < total_length:
			if seqname + ":" + str(pos) not in avoid:
				outfile.write(seqname + "\t" + str(pos) + "\t" + str(vaf)+ "\n")
				vaf += vafspacing
				if vaf >= 1:
					vaf = vafspacing
			pos += variantspacing
		interval = faifile.readline().strip().split()
