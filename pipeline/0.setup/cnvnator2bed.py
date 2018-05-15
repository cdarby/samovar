import sys
#python cnvnator2bed.py [infile] [outfile] [.fai of genome]

if len(sys.argv) < 4: sys.exit()
faiDict = dict()
with open(sys.argv[3]) as fai:
	L = fai.readline().split()
	while len(L) > 0:
		try:
			faiDict[L[0]] = int(L[1])
			L = fai.readline().split()
		except ValueError:
			L = fai.readline().split()

with open(sys.argv[1]) as F, open(sys.argv[2],"w") as out:
	L = F.readline().split()
	#assume that cnvnator coordinates are 1-based like vcf file -> bed is 0-based
	while len(L) > 0:
		q = L[1]
		q = q.replace(":"," ")
		q = q.replace("-"," ")
		q = q.split()
		if len(q) < 2: print(q)
		if int(q[1]) < 6: 
			lo = 0
		else:
			lo = int(q[1]) - 6
		hi = int(q[2]) + 5
		if hi > faiDict[q[0]]: hi = faiDict[q[0]]
		out.write(q[0] + "\t" + str(lo) + "\t" + str(hi) + "\n")

		L = F.readline().split()