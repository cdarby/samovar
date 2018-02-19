
import pysam,numpy,argparse

NT = "ACGTN"

parser = argparse.ArgumentParser(description='python siteFeatures.py --bam diploid_output/phased_possorted_bam.bam --bed haploidPos.bed --out diploidBam_haploidPos.pileup')
parser.add_argument('--bam', help='Input bam file name',required=True)
parser.add_argument('--bed', help='Input bed file name',required=True)
parser.add_argument('--out', help='Output file name',required=True)
args = parser.parse_args()


R = []
#Regions
with open(args.bed) as bedfile:
	region = bedfile.readline().strip().split()
	while len(region) > 0:
		R.append(region)
		region = bedfile.readline().strip().split()

outfile = open(args.out + args.bed +".filteringfeatures", "w")
samfile = pysam.AlignmentFile(args.bam, "rb")

for i in xrange(len(R)):
	region = R[i]
	chrom = region[0]
	start = int(region[1])# -1 # -1 for MosaicHunter
	end = start + 1 
	#print region
	for pileupcolumn in samfile.pileup(chrom, start, end): 
		#allreads which overlap the region are returned. The first base returned will be the first base of the first read not necessarily the first base of the region used in the query.
		if pileupcolumn.pos < start:
			continue
		if pileupcolumn.pos >= end:
			break

		readpos = []
		bases = []
		clip = []
		ind = []

		XH = []
		XL = []
		XV = []
		XR = []
		XB = []
		
		depth = pileupcolumn.nsegments
		for pileupread in pileupcolumn.pileups:
			aln = pileupread.alignment

			#base and base quality
			if pileupread.is_del and not pileupread.is_refskip:
				depth -= 1
				continue     
			else:
				bases.append(str(aln.query_sequence[pileupread.query_position]))
			

			#position on read
			#readpos.append(pileupread.query_position_or_next) #result of this is strange
			readpos.append(min((aln.query_alignment_end - pileupread.query_position_or_next),(pileupread.query_position_or_next - aln.query_alignment_start)))

			#fraction of alignment positions that are soft+hardclip
			#g  = aln.get_cigar_stats()[0]
			#clip.append(1.0*(g[4]+g[5])/sum(g[:10]),3)
			#query_alignment_sequence: This is a substring of seq that excludes flanking bases that were soft clipped 
			clip.append(aln.query_length - len(aln.query_alignment_sequence))


			#fraction of alignment positions that are indel
			#ind.append(1.0*(g[1]+g[2])/sum(g[:10]),3)
			#ind.append(abs(aln.query_length-aln.reference_length))
			n_ind = 0
			for (op,ln) in aln.cigartuples:
				if op == 1 or op == 2 or op == 3:
					n_ind += ln
			ind.append(n_ind)

			# linked read-quality features
			if aln.has_tag("HP"): #Longranger read phasing tag
				XH.append(str(aln.get_tag("HP")-1))
			else: 
				XH.append("*")
			if aln.has_tag("XL"):
				XL.append(round(aln.get_tag("XL"),3))
			else:
				XL.append("*")
			#if aln.has_tag("XZ"):
			#    XZ.append(round(aln.get_tag("XZ"),3))
			#else:
			#    XZ.append(0.0) 
			if aln.has_tag("XV"):
				XV.append(aln.get_tag("XV"))
				XR.append(aln.get_tag("XR"))
				XB.append(aln.get_tag("XB"))
			else:
				XV.append(0)
				XR.append(0) 
				XB.append(0)

		if depth <= 0: 
			continue

		#features for whole site
		nphased = sum([1 for i in xrange(depth) if XH[i] != "*"])
		fracphased = str(round(1.0*nphased/depth,3))
		basecounts = [bases.count(B) for B in NT]
		basecount_indexes = [i[0] for i in sorted(enumerate(basecounts), key=lambda x:x[1])]
		major = NT[basecount_indexes[-1]] #most frequent nt
		minor = NT[basecount_indexes[-2]]
		if minor == "N" : minor == NT[basecount_indexes[-3]] #doesn't seem to be working
		#remove all reads that are not major or minor
		depth_orig = depth
		oldlists = [bases,readpos,XH,XL,XV,XR,XB,clip,ind]
		newlists = []
		for L in oldlists:
			newlists.append([L[i] for i in xrange(depth) if bases[i] == minor or bases[i] == major])
		depth = len(newlists[0])
		(bases,readpos,XH,XL,XV,XR,XB,clip,ind) = newlists
		fracOtherAllele = str(round(1.0*(depth_orig-depth)/depth,3))
		nminor = sum([1 for i in xrange(depth) if bases[i] != major])
		if nminor == 0: continue

		frach0 = 0.0 if nphased == 0 else 1.0*sum([1 for i in xrange(depth) if XH[i] == "0"])/nphased
		frach1 = 0.0 if nphased == 0 else 1.0*sum([1 for i in xrange(depth) if XH[i] == "1"])/nphased
		MAF_h0 = "0" if nphased == 0 else str(round(1.0*sum([1 for i in xrange(depth) if XH[i] == "0" and bases[i] == minor])/nphased ,3))
		#print 1.0*sum([1 for i in xrange(depth) if XH[i] == "0" and bases[i] == minor])/nphased
		MAF_h1 = "0" if nphased == 0 else str(round(1.0*sum([1 for i in xrange(depth) if XH[i] == "1" and bases[i] == minor])/nphased ,3))

		frach = str(round(max(frach0,frach1),3))
		MAF = str(round(1.0*sum([1 for i in xrange(depth) if bases[i] != major])/depth,3)) #of all reads
		MAF_phased = str(0.0 if nphased == 0 else round(1.0*sum([1 for i in xrange(depth) if bases[i] != major and XH[i] != "*"])/nphased,3)) #of phased reads

		outfile.write("\t".join([chrom,str(start),str(depth),fracphased,frach,MAF,MAF_phased,MAF_h0,MAF_h1]))
		outfile.write("\t")

		#linked features

		#Reads that need to be changed in the allele x haplotype matrix
		change_to_hom = sum([1 for i in xrange(depth) if XH[i] != "*" and bases[i] != major])
		change_to_het_left = sum([1 for i in xrange(depth) if (XH[i] == "1" and bases[i] == major) or (XH[i] == "0" and bases[i] != major)])
		change_to_het_right = sum([1 for i in xrange(depth) if (XH[i] == "0" and bases[i] == major) or (XH[i] == "1" and bases[i] != major)])
		
		#ones to Change
		oldlists = [bases,readpos,XH,XL,XV,XR,XB,clip,ind]
		Changelists = []
		Notchangelists = []

		if change_to_hom <= change_to_het_left and change_to_hom <= change_to_het_right:
			for L in oldlists:
				Changelists.append([L[i] for i in xrange(depth) if XH[i] != "*" and bases[i] != major])
				Notchangelists.append([L[i] for i in xrange(depth) if XH[i] != "*" and bases[i] == major])
		elif change_to_het_left <= change_to_hom  and change_to_het_left <= change_to_het_right:
			#print "change_to_het_left"
			for k,L in enumerate(oldlists):
				Changelists.append([L[i] for i in xrange(depth) if (XH[i] == "1" and bases[i] == major) or (XH[i] == "0" and bases[i] != major)])
				Notchangelists.append([L[i] for i in xrange(depth) if (XH[i] == "0" and bases[i] == major) or (XH[i] == "1" and bases[i] != major)])
		else:   
			for k,L in enumerate(oldlists):
				#print "change_to_het_right"
				Changelists.append([L[i] for i in xrange(depth) if (XH[i] == "0" and bases[i] == major) or (XH[i] == "1" and bases[i] != major)])
				Notchangelists.append([L[i] for i in xrange(depth) if (XH[i] == "1" and bases[i] == major) or (XH[i] == "0" and bases[i] != major)])

		(Nbases,Nreadpos,NXH,NXL,NXV,NXR,NXB,Nclip,Nind) = Notchangelists
		(Cbases,Creadpos,CXH,CXL,CXV,CXR,CXB,Cclip,Cind) = Changelists
		if len(Cbases) == 0:
			outfile.write("\t".join(["0"]*12))
			outfile.write("\t")
			nC = 0
		else:
			nC = len(Cbases)
			Cfrach0= 1.0*sum([1 for i in xrange(nC) if CXH[i] == "0"])/nC
			Cfrach1= 1.0*sum([1 for i in xrange(nC) if CXH[i] == "1"])/nC
			Cfrach= str(round(max(Cfrach0,Cfrach1),3))
			Cbasecounts = [Cbases.count(B) for B in NT]
			Cmajor = NT[Cbasecounts.index(max(Cbasecounts))] #most frequent nt
			CMAF = str(round(1.0*sum([1 for i in xrange(nC) if Cbases[i] != Cmajor])/nC,3)) #of all reads
			CavgXL = "0" #str(round(1.0*sum(CXL) / nC,3))
			CavgXV = "0" #str(round(1.0*sum(CXV) / nC,3))
			CavgXR = "0" #str(round(1.0*sum(CXR) / nC,3))
			CavgXB = "0" #str(round(1.0*sum(CXB) / nC,3))
			Cavgpos = str((round(1.0*sum(Creadpos) / nC,3)))
			Cavgclip = str((round(1.0*sum(Cclip) / nC,3)))
			Cavgind = str((round(1.0*sum(Cind) / nC,3)))
			Cfracclip = str(round(1.0*sum([1 if i>0 else 0 for i in Cclip]) / nC,3))
			Cfracind = str(round(1.0*sum([1 if i>0 else 0 for i in Cind]) / nC,3))
			fracC = str(round(1.0*nC/(nphased),3))

			outfile.write("\t".join([fracC,Cfrach,CMAF,CavgXL,CavgXV,CavgXR,CavgXB,Cavgpos,Cavgclip,Cfracclip,Cavgind,Cfracind]))
			outfile.write("\t")

		if len(Nbases) == 0:
			outfile.write("\t".join(["0"]*10))
			outfile.write("\t")
			nN = 0
		else:
			nN = len(Nbases)
			Nfrach0= 1.0*sum([1 for i in xrange(nN) if NXH[i] == "0"])/nN
			Nfrach1= 1.0*sum([1 for i in xrange(nN) if NXH[i] == "1"])/nN
			Nfrach= str(round(max(Nfrach0,Nfrach1),3))
			Nbasecounts = [Nbases.count(B) for B in NT]
			Nmajor = NT[Nbasecounts.index(max(Nbasecounts))] #most frequent nt
			NMAF = str(round(1.0*sum([1 for i in xrange(nN) if Nbases[i] != Nmajor])/nN,3)) #of all reads
			NavgXL = "0" #str(round(1.0*sum(NXL) / nN,3))
			NavgXV = "0" #str(round(1.0*sum(NXV) / nN,3))
			NavgXR = '0' #str(round(1.0*sum(NXR) / nN,3))
			NavgXB = '0' #str(round(1.0*sum(NXB) / nN,3))
			Navgpos = str(round(1.0*sum(Nreadpos) / nN,3))
			Navgclip = str((round(1.0*sum(Nclip) / nN,3)))
			Navgind = str((round(1.0*sum(Nind) / nN,3)))
			Nfracclip = str(round(1.0*sum([1 if i>0 else 0 for i in Nclip]) / nN,3))
			Nfracind = str(round(1.0*sum([1 if i>0 else 0 for i in Nind]) / nN,3))

			outfile.write("\t".join([Nfrach,NMAF,NavgXL,NavgXV,NavgXR,NavgXB,Navgpos,Navgclip,Nfracclip,Navgind,Nfracind]))
			outfile.write("\t")

		#non-linked features - "M-reads" are those with minor allele (like C-reads), "J-reads" are those with major allele (like N-reads)
		nminor = sum([1 for i in xrange(depth) if bases[i] != major])
		Mavgpos = str(round(1.0*sum([readpos[i] for i in xrange(depth) if bases[i] != major]) / nminor,3))
		Mavgclip = str(round(1.0*sum([clip[i] for i in xrange(depth) if bases[i] != major]) / nminor,3))
		Mavgind = str(round(1.0*sum([ind[i] for i in xrange(depth) if bases[i] != major]) / nminor,3))
		Mclip = [clip[i] for i in xrange(depth) if bases[i] != major]
		Mfracclip = str(round(1.0*sum([1 if i>0 else 0 for i in Mclip]) / nminor,3))
		Mind = [ind[i] for i in xrange(depth) if bases[i] != major]
		Mfracind = str(round(1.0*sum([1 if i>0 else 0 for i in Mind]) / nminor,3))

		nmajor = depth-nminor
		Javgpos = str(round(1.0*sum([readpos[i] for i in xrange(depth) if bases[i] == major]) / nmajor,3))
		Javgclip = str(round(1.0*sum([clip[i] for i in xrange(depth) if bases[i] == major]) / nmajor,3))
		Javgind = str(round(1.0*sum([ind[i] for i in xrange(depth) if bases[i] == major]) / nmajor,3))
		
		Jind = [ind[i] for i in xrange(depth) if bases[i] == major]
		Jfracind = str(round(1.0*sum([1 if i>0 else 0 for i in Jind]) / nmajor,3))		
		Jclip = [clip[i] for i in xrange(depth) if bases[i] == major]
		Jfracclip = str(round(1.0*sum([1 if i>0 else 0 for i in Jclip]) / nmajor,3))
		MAFnorm = str(round(abs(0.25-1.0*nminor/depth),3))

		outfile.write("\t".join([Mavgpos,Mavgclip,Mfracclip,Mavgind,Mfracind,Javgpos,Javgclip,Jfracclip,Javgind,Jfracind,MAFnorm]))
		outfile.write("\n")


samfile.close()
outfile.close()
