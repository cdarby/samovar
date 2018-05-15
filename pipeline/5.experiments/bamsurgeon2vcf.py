import sys
# bamsurgeon2vcf.py < [infile]

L = sys.stdin.readline().strip().split()
while len(L) > 0:
	chr = L[1].split(":")[0]
	pos = L[3]
	ref = L[4][0]
	alt = L[4][-1]
	print(chr + "\t" + pos + "\t.\t" + ref + "\t" + alt + "\t0\tPASS\t.\tGT\t0|0")
	L = sys.stdin.readline().strip().split()