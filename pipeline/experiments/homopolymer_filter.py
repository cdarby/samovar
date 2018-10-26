import pysam

F = pysam.FastaFile("/work-zfs/mschatz1/resources/refdata-hg19-2.1.0/fasta/genome.fa")
N_WINDOWS = 0
NUM_BASES_MASKED = 0
with open("~/scratch/MosaicHunter/resources/hg38.complement.bed") as bedfile:
	L = bedfile.readline().strip().split()
	while len(L) > 0:
		N_WINDOWS += 1
		seq = F.fetch(reference=region[0],start=int(region[1]),end=int(region[2])+1)
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
					space_btw = (i - hp_len + 1) - last_hp_end + 1
					NUM_BASES_MASKED += hp_len + min(space_btw - 2, 2) + 2
					last_hp_end = i
				elif hp_len > 5:
					space_btw = (i - hp_len + 1) - last_hp_end + 1
					NUM_BASES_MASKED += hp_len + min(space_btw - 3, 3) + 3
					last_hp_end = i
				hp_len = 1
			#if long hp have less than 6bp between or short hp have less than 4bp between, not considering if there is at least that much space between overcounts the number of bases masked
		if N_WINDOWS % 50000 == 0: print(NUM_BASES_MASKED)
		L = bedfile.readline().strip().split()

print(NUM_BASES_MASKED)

#MH Visible: complement of intersect of WGS regions; all_repeats; homopolymer
#b37: REMAINING 1469041925 - 267822039 masked = 1201219886 visible ~39%
#hg38: REMAINING 1497984122 - 392797495 masked = 1105186627 visible

#Samovar Visible: complement of repeatasm bedfile
#b37: LOST 935652587 - TOTAL 3095677412 = 2160024825 visible 69.8
#hg38: LOST 878411944 - TOTAL 3088269832 = 2209857888 visible 71.6

#HM Visible: copmlement of self alignments; STR [missing BLAT]

#All Visible (except BLAT)
#b37: 725095971
#hg38: 695174817






