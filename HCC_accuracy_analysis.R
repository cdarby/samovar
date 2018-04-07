

#### Only repeat region filter - no feature based filters

colnamesShort <- c("status","id","chrom","pos","score","depth","fracphased","fracstrand","frach","MAF","MAF_phased",
              "nC","fracC","Cavgbq","Cfrach","Cavgmapq","Cfracstrand","CMAF","CavgXL","CavgXV","CavgXR","CavgXB","Cavgpos","Cavgclip","Cavgind","Cfracclip","Cfracind","CavgASXS","Cavgpair",
              "Navgbq","Nfrach","Navgmapq","Nfracstrand","NMAF","NavgXL","NavgXV","NavgXR","NavgXB","Navgpos","Navgclip","Nfracclip","Nfracind","Navgind","NavgASXS","Navgpair",
              "weightedCbq","Cavgibq",
              "Mavgbq","Mavgmapq","Mfracstrand","Mavgpos","Mavgclip","Mavgind","Mfracclip","Mfracind","Mavgibq","weightedMbq","MavgASXS","Mavgpair",
              "Javgbq","Javgmapq","Jfracstrand","Javgpos","Javgclip","Javgind","Jfracclip","Jfracind","JavgASXS","Javgpair","MAFnorm","fracOtherAllele")

colnamesLong <- c("status","id","chromlinkage","poslinkage","scorelinkage","scorelinkagead","linkagemat","chrom","pos","score","depth","fracphased","fracstrand","frach","MAF","MAF_phased",
              "nC","fracC","Cavgbq","Cfrach","Cavgmapq","Cfracstrand","CMAF","CavgXL","CavgXV","CavgXR","CavgXB","Cavgpos","Cavgclip","Cavgind","Cfracclip","Cfracind","CavgASXS","Cavgpair",
              "Navgbq","Nfrach","Navgmapq","Nfracstrand","NMAF","NavgXL","NavgXV","NavgXR","NavgXB","Navgpos","Navgclip","Nfracclip","Nfracind","Navgind","NavgASXS","Navgpair",
              "weightedCbq","Cavgibq",
              "Mavgbq","Mavgmapq","Mfracstrand","Mavgpos","Mavgclip","Mavgind","Mfracclip","Mfracind","Mavgibq","weightedMbq","MavgASXS","Mavgpair",
              "Javgbq","Javgmapq","Jfracstrand","Javgpos","Javgclip","Javgind","Jfracclip","Jfracind","JavgASXS","Javgpair","MAFnorm","fracOtherAllele")


NB <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/completefeatures/scanNB_norepeat_labeled.features.tsv",sep="\t",header=F,col.names=colnamesShort)
NE <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/completefeatures/scanNE_norepeat_labeled.features.tsv",sep="\t",header=F,col.names=colnamesShort)
TB <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/completefeatures/scanTB_norepeat_linkage_labeled.features.tsv",sep="\t",header=F,col.names=colnamesLong)
TE <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/completefeatures/scanTE_norepeat_linkage_labeled.features.tsv",sep="\t",header=F,col.names=colnamesLong)

# NB <- NB[order(NB$V5,decreasing=T),]
# NE <- NE[order(NE$V5,decreasing=T),]
# TB <- TB[order(TB$V5,decreasing=T),]
# TE <- TE[order(TE$V5,decreasing=T),]

#real status = 2; sim status = 1
TE_noB <- TE[TE$status != 1,]
TB_noB <- TB[TB$status != 1,]

TE_noR <- TE[TE$status != 2,]
TB_noR <- TB[TB$status != 2,]

# length(TE$status[TE$status == 0])
# length(TE$status[TE$status == 1])
# length(TE$status[TE$status == 2])
# hist(TE$score[TE$status == 0])
# hist(TE$score[TE$status == 1])
# hist(TE$score[TE$status == 2])
# hist(TE$MAF[TE$status == 0])
# hist(TE$MAF[TE$status == 1])
# hist(TE$MAF[TE$status == 2])
# 
# length(TE_noB$status[TE_noB$status == 0])
# length(TE_noB$status[TE_noB$status == 2])
# hist(TE_noB$score[TE_noB$status == 0])
# hist(TE_noB$score[TE_noB$status == 2])
# hist(TE_noB$MAF[TE_noB$status == 0])
# hist(TE_noB$MAF[TE_noB$status == 2])

par(mfrow=c(4,3))
for (i in seq(0,0.45,by=0.05)) {
  t <- TB_noR[TB_noR$MAF >= i & TB_noR$MAF < i+0.05 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  plot(main=paste(i,length(t$status), round(auc(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction)),5)),getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
}


calls <- TE_noB
filteredTE <- calls[calls$scorelinkage > 0.05 & calls$fracC>0 & calls$fracphased>0 & calls$Cfracclip<1 & calls$Mfracclip<1 & calls$Nfracclip<1 
                    & calls$Jfracclip<1 & calls$frach != 0 & calls$frach != 1 & calls$Cfracind<1 & calls$Mfracind<1 & calls$Nfracind<1 & calls$Jfracind<1 
                    #& (calls$Cfrach > 0.9 | calls$Cfrach < 0.1) 
                    & calls$Mavgpos>10 & calls$Cavgpos>10 & calls$Javgpos>10 
                    #& calls$frach<0.7 
                    & calls$frach > 0.3 & calls$Cavgpair < 1 & calls$Mavgpair < 1 & calls$Navgpair < 1 & calls$Javgpair < 1 &  calls$Cfracstrand < 1 
                    & calls$Cfracstrand > 0 & calls$Mfracstrand < 1 & calls$Mfracstrand > 0 & calls$Nfracstrand < 1 & calls$Nfracstrand > 0 
                    & calls$Jfracstrand < 1 & calls$Jfracstrand > 0 & calls$MAF>0.05 
                    #& calls$fracphased>0.7 
                    & calls$fracOtherAllele < 0.05 
                    & calls$depth>=16 & calls$nC>=4,]
calls <- TB_noB
filteredTB <- calls[calls$scorelinkage > 0.05 & calls$fracC>0 & calls$fracphased>0 & calls$Cfracclip<1 & calls$Mfracclip<1 & calls$Nfracclip<1 
                    & calls$Jfracclip<1 & calls$frach != 0 & calls$frach != 1 & calls$Cfracind<1 & calls$Mfracind<1 & calls$Nfracind<1 & calls$Jfracind<1 
                    #& (calls$Cfrach > 0.9 | calls$Cfrach < 0.1) 
                    & calls$Mavgpos>10 & calls$Cavgpos>10 & calls$Javgpos>10 
                    #& calls$frach<0.7 
                    & calls$frach > 0.3 & calls$Cavgpair < 1 & calls$Mavgpair < 1 & calls$Navgpair < 1 & calls$Javgpair < 1 &  calls$Cfracstrand < 1 
                    & calls$Cfracstrand > 0 & calls$Mfracstrand < 1 & calls$Mfracstrand > 0 & calls$Nfracstrand < 1 & calls$Nfracstrand > 0 
                    & calls$Jfracstrand < 1 & calls$Jfracstrand > 0 & calls$MAF>0.05 
                    #& calls$fracphased>0.7 
                    & calls$fracOtherAllele < 0.05 
                    & calls$depth>=16 & calls$nC>=4,]


par(mfrow=c(3,3))
for (i in seq(0.05,0.45,by=0.05)) {
  t <- filteredTB[filteredTB$MAF >= i & filteredTB$MAF < i+0.05 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  plot(main=paste(i,length(t$status),round(auc(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction)),5)),getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
}

#auc(getaxis(filteredTB,A,xaxisfunction),getaxis(filteredTB,A,yaxisfunction))

K <- colorRampPalette(c("red","purple","cyan"))
m <- K(9)
S <- seq(0.05,0.45,by=0.05)

plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="tumor/edited-training (bamsurgeon sites; linkage filter ONLY)")
c <- 1
for (i in S) {
  t <- TE_noR[TE_noR$MAF >= i & TE_noR$MAF < i+0.05 & TE_noR$scorelinkage > 0.05,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
}
B <- seq(min(TE_noR$score),max(TE_noR$score),length.out = 100)
lines(getaxis(TE_noR,B,xaxisfunction),getaxis(TE_noR,B,yaxisfunction))
legend("bottomright",legend=S,col=m,lty=rep(1,9)) 


#### Repeat region and feature based (strict?) filters applied

colnames <- c("chromlinkage","poslinkage","scorelinkage","scorelinkagead","linkagemat","chrom","pos","score","depth","fracphased","fracstrand","frach","MAF","MAF_phased",
              "nC","fracC","Cavgbq","Cfrach","Cavgmapq","Cfracstrand","CMAF","CavgXL","CavgXV","CavgXR","CavgXB","Cavgpos","Cavgclip","Cavgind","Cfracclip","Cfracind","CavgASXS","Cavgpair",
              "Navgbq","Nfrach","Navgmapq","Nfracstrand","NMAF","NavgXL","NavgXV","NavgXR","NavgXB","Navgpos","Navgclip","Nfracclip","Nfracind","Navgind","NavgASXS","Navgpair",
              "weightedCbq","Cavgibq",
              "Mavgbq","Mavgmapq","Mfracstrand","Mavgpos","Mavgclip","Mavgind","Mfracclip","Mfracind","Mavgibq","weightedMbq","MavgASXS","Mavgpair",
              "Javgbq","Javgmapq","Jfracstrand","Javgpos","Javgclip","Javgind","Jfracclip","Jfracind","JavgASXS","Javgpair","MAFnorm","fracOtherAllele")


NB <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCNBS.linkage.tsv",sep="\t",header=F,col.names=colnames)
NE <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCNE.linkage.tsv",sep="\t",header=F,col.names=colnames)
TB <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTBS.linkage.tsv",sep="\t",header=F,col.names=colnames)
TE <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTE.linkage.tsv",sep="\t",header=F,col.names=colnames)


length(NB$V3[NB$V3>0.05])
# 3567
length(NE$V3[NE$V3>0.05])
# 3450
length(TE$V3[TE$V3>0.05])
# 1675
length(TB$V3[TB$V3>0.05])
# 1686

hist(NB$V13)
hist(NE$V13)
hist(TE$V13)
hist(TB$V13)

NV <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/sitesNvisible.linkage.tsv",sep="\t",header=F)
TV <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/sitesTvisible.linkage.tsv",sep="\t",header=F)
length(NV$V3[NV$V3>0.05])
#
length(TV$V3[TV$V3>0.05])
#1105

write.table(NB[NB$V3 > 0.05,],file="~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCNBS.linkagefiltered.tsv",row.names=F,col.names=F,sep="\t",quote=F) 
write.table(NE[NE$V3 > 0.05,],file="~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCNE.linkagefiltered.tsv",row.names=F,col.names=F,sep="\t",quote=F) 
write.table(TB[TB$V3 > 0.05,],file="~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTBS.linkagefiltered.tsv",row.names=F,col.names=F,sep="\t",quote=F) 
write.table(TE[TE$V3 > 0.05,],file="~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTE.linkagefiltered.tsv",row.names=F,col.names=F,sep="\t",quote=F) 


### Old (more strict) feature based filters were applied 

A<-seq(0.9,1,0.001)
par(mfrow=c(2,2))
TE <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTE_labeled.txt",sep="\t",header=F,col.names=c("status","pos","score"))
TB <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTBS_labeled.txt",sep="\t",header=F,col.names=c("status","pos","score"))
NE <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCNE_labeled.txt",sep="\t",header=F,col.names=c("status","pos","score"))
NB <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCNBS_labeled.txt",sep="\t",header=F,col.names=c("status","pos","score"))

plot(main=paste("HCCT-editedsites",round(auc(getaxis(TE,A,xaxisfunction),getaxis(TE,A,yaxisfunction)),5)),getaxis(TE,A,xaxisfunction),getaxis(TE,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
plot(main=paste("HCCT-bamsurgeon",round(auc(getaxis(TB,A,xaxisfunction),getaxis(TB,A,yaxisfunction)),5)),getaxis(TB,A,xaxisfunction),getaxis(TB,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
plot(main=paste("HCCN-editedsites",round(auc(getaxis(NE,A,xaxisfunction),getaxis(NE,A,yaxisfunction)),5)),getaxis(NE,A,xaxisfunction),getaxis(NE,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
plot(main=paste("HCCN-bamsurgeon",round(auc(getaxis(NB,A,xaxisfunction),getaxis(NB,A,yaxisfunction)),5)),getaxis(NB,A,xaxisfunction),getaxis(NB,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))

#real status = 1; sim status = 9
TE1 <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTE_labeled1.txt",sep="\t",header=F,col.names=c("status","pos","score"))
TE1 <- TE1[TE1$status < 2,]
plot(main=paste("HCCT-editedsites (real)",round(auc(getaxis(TE1,A,xaxisfunction),getaxis(TE1,A,yaxisfunction)),5)),getaxis(TE1,A,xaxisfunction),getaxis(TE1,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))

TE1 <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTE_labeled1.txt",sep="\t",header=F,col.names=c("status","pos","score"))
TE1 <- TE1[TE1$status != 1,]
plot(main=paste("HCCT-editedsites (sim)",round(auc(getaxis(TE1,A,xaxisfunction),getaxis(TE1,A,yaxisfunction)),5)),getaxis(TE1,A,xaxisfunction),getaxis(TE1,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))

TB1 <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTBS_labeled1.txt",sep="\t",header=F,col.names=c("status","pos","score"))
TB1 <- TB1[TB1$status < 2,]
plot(main=paste("HCCT-bamsurgeon (real)",round(auc(getaxis(TB1,A,xaxisfunction),getaxis(TB1,A,yaxisfunction)),5)),getaxis(TB1,A,xaxisfunction),getaxis(TB1,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))

TB1 <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/HCCTBS_labeled1.txt",sep="\t",header=F,col.names=c("status","pos","score"))
TB1 <- TB1[TB1$status != 1,]
plot(main=paste("HCCT-bamsurgeon (sim)",round(auc(getaxis(TB1,A,xaxisfunction),getaxis(TB1,A,yaxisfunction)),5)),getaxis(TB1,A,xaxisfunction),getaxis(TB1,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))


###################

rank <- seq(1,length(filteredTE$status))
S <- seq(0,0.45,by=0.05)
filteredTE <- filteredTE[order(filteredTE$score,decreasing=T),]
sum(filteredTE$status[1:1000])/2

sapply(S,function(i) length(filteredTE$status[filteredTE$MAF >= i & filteredTE$MAF < i+0.05]))
sapply(S,function(i) sum(filteredTE$status[filteredTE$MAF >= i & filteredTE$MAF < i+0.05])/2)

sapply(S,function(i) sum(filteredTE$status[filteredTE$MAF >= i & filteredTE$MAF < i+0.05][1:100])/2)
sapply(S,function(i) length(filteredTE$MAF[1:10][filteredTE$MAF[1:10] >= i & filteredTE$MAF[1:10] < i+0.05]))
sapply(S,function(i) mean(filteredTE$score[filteredTE$MAF >= i & filteredTE$MAF < i+0.05][1:10]))
sapply(S,function(i) mean(rank[filteredTE$MAF >= i & filteredTE$MAF < i+0.05][1:10]))








