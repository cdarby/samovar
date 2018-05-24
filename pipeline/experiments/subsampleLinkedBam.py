import pysam,random,sys
#python subsampleLinkedBam.py [infile] [outfile]
if len(sys.argv) < 3: sys.exit()

samfile = pysam.AlignmentFile(sys.argv[2], "rb")
outfile = pysam.AlignmentFile(sys.argv[3],"wb",template=samfile)
for read in samfile.fetch():
	if read.has_tag("BX"):
		t = read.get_tag("BX")
		if t[9] == 'A' or t[9] == 'G':
			outfile.write(read)
	else:
		if random.random() > 0.5:
			outfile.write(read)
samfile.close()
outfile.close()