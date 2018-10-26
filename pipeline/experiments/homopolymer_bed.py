import pysam,sys

F = pysam.FastaFile("/work-zfs/mschatz1/resources/refdata-GRCh38-2.1.0/fasta/genome.fa")

with open("hg38.txt") as bedfile:
	L = bedfile.readline().strip().split()
	while len(L) > 0:
		s = int(L[1])
		e = int(L[2])+1
		c = L[0]
		seq = F.fetch(reference=c,start=s,end=e)
		#A region can either be specified by reference, start and end. start and end denote 0-based, half-open intervals.
		hp_char = seq[0]
		hp_len = 1
		last_hp_end = 0

		for i in range(1,len(seq)):
			if seq[i] == hp_char:
				hp_len += 1
			else:
				hp_char = seq[i]
				if hp_len == 4 or hp_len == 5:
					print(c + "\t" + str(min(0,s+i-hp_len-2)) + "\t" + str(max(e,s+i+2))) #buffer of 2
				elif hp_len > 5:
					print(c + "\t" + str(min(0,s+i-hp_len-3)) + "\t" + str(max(e,s+i+3))) #buffer of 3
				hp_len = 1
		L = bedfile.readline().strip().split()

'''
with open("selfchain_hg19_main.bed") as bedfile:
	L = bedfile.readline().strip().split()
	while len(L) > 0:
		if int(L[1]) < int(L[2]):
			print("\t".join(L))
		else:
			print("\t".join([L[0], L[2], L[1]]))
			L = bedfile.readline().strip().split()
'''


