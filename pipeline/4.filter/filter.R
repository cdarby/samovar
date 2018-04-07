#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "featurefiltered.tsv"
}

colnames <- c("chrom","pos","score","depth","fracphased","fracstrand","frach","MAF","MAF_phased",
              "nC","fracC","Cavgbq","Cfrach","Cavgmapq","Cfracstrand","CMAF","CavgXL","CavgXV","CavgXR","CavgXB","Cavgpos","Cavgclip","Cavgind","Cfracclip","Cfracind","CavgASXS","Cavgpair",
              "Navgbq","Nfrach","Navgmapq","Nfracstrand","NMAF","NavgXL","NavgXV","NavgXR","NavgXB","Navgpos","Navgclip","Nfracclip","Nfracind","Navgind","NavgASXS","Navgpair",
              "weightedCbq","Cavgibq",
              "Mavgbq","Mavgmapq","Mfracstrand","Mavgpos","Mavgclip","Mavgind","Mfracclip","Mfracind","Mavgibq","weightedMbq","MavgASXS","Mavgpair",
              "Javgbq","Javgmapq","Jfracstrand","Javgpos","Javgclip","Javgind","Jfracclip","Jfracind","JavgASXS","Javgpair","MAFnorm","fracOtherAllele")

calls <- read.delim(args[1],header=FALSE,sep="\t",col.names=colnames)

filtered <- calls[calls$fracC>0 & calls$fracphased>0 & calls$Cfracclip<1 & calls$Mfracclip<1 & calls$Nfracclip<1 & calls$Jfracclip<1 & calls$frach != 0 & calls$frach != 1 & calls$Cfracind<1 & calls$Mfracind<1 & calls$Nfracind<1 & calls$Jfracind<1 & (calls$Cfrach > 0.9 | calls$Cfrach < 0.1) & calls$Mavgpos>10 & calls$Cavgpos>10 & calls$Javgpos>10 & calls$frach<0.7 & calls$frach > 0.3 & calls$Cavgpair < 1 & calls$Mavgpair < 1 & calls$Navgpair < 1 & calls$Javgpair < 1 &  calls$Cfracstrand < 1 & calls$Cfracstrand > 0 & calls$Mfracstrand < 1 & calls$Mfracstrand > 0 & calls$Nfracstrand < 1 & calls$Nfracstrand > 0 & calls$Jfracstrand < 1 & calls$Jfracstrand > 0 & calls$MAF>0.05 & calls$fracphased>0.7 & calls$fracOtherAllele < 0.05 & calls$depth>=16 & calls$nC>=4,]

write.table(filtered,file=args[2],row.names=F,col.names=F,sep="\t",quote=F) 
#& calls$Navgpos>10 