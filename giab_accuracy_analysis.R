###### Paired, region filter


colnamesPaired <- c("status","id","pairedscore","chrom","pos","unpairedscore","depth","MAF","Mavgbq","Mavgpos","Mavgclip","Mavgind","weightedMbq","Javgbq","Javgpos","Javgclip","Javgind")
                    #,"fracphased","frach","MAF_phased","nC","fracC","Cavgbq","Cfrach","Cavgmapq","Cfracstrand","CMAF","Cavgpos","Cavgclip","Cavgind","Cfracclip","Cfracind","Cavgpair",
                    #"Navgbq","Nfrach","Navgmapq","Nfracstrand","NMAF","NavgXL","NavgXV","NavgXR","NavgXB","Navgpos","Navgclip","Nfracclip","Nfracind","Navgind","NavgASXS","Navgpair","weightedCbq")
pairedend <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/pairedend/paired.features_labeled.tsv",sep="\t",header=F,col.names=colnamesPaired)

pairedend$score <- pairedend$unpairedscore
K <- colorRampPalette(c("red","purple","cyan"))
m <- K(10)
S <- seq(0,0.45,by=0.05)

plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="No Phasing - repeat filter only")
c <- 1
for (i in S) {
  t <- pairedend[pairedend$MAF >= i & pairedend$MAF < i+0.05 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
} 
B <- seq(min(pairedend$score),max(pairedend$score),length.out = 100)
lines(getaxis(pairedend,B,xaxisfunction),getaxis(pairedend,B,yaxisfunction),type="l")
legend("bottomright",legend=S,col=m,lty=rep(1,10)) 



pairedend$score <- pairedend$pairedscore
plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="Only sites with paired end data - repeat filter only")
c <- 1
for (i in S) {
  t <- pairedend[pairedend$MAF >= i & pairedend$MAF < i+0.05 & pairedend$score > 0 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
} 
B <- seq(min(pairedend$score),max(pairedend$score),length.out = 100)
lines(getaxis(pairedend[pairedend$score > 0,],B,xaxisfunction),getaxis(pairedend[pairedend$score > 0,],B,yaxisfunction),type="l")
legend("bottomright",legend=S,col=m,lty=rep(1,10)) 



pairedend$score <- pairedend$pairedscore
plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="Only sites with paired end data - repeat filter only")
c <- 1
for (i in S) {
  t <- pairedend[pairedend$MAF >= i & pairedend$MAF < i+0.05 & pairedend$score > 0 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
} 
B <- seq(min(pairedend$score),max(pairedend$score),length.out = 100)
lines(getaxis(pairedend[pairedend$score > 0,],B,xaxisfunction),getaxis(pairedend[pairedend$score > 0,],B,yaxisfunction),type="l")
legend("bottomright",legend=S,col=m,lty=rep(1,10)) 


pairedend$score <- pmax(pairedend$pairedscore,pairedend$unpairedscore)
plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="Hybrid scoring - repeat filter only")
c <- 1
for (i in S) {
  t <- pairedend[pairedend$MAF >= i & pairedend$MAF < i+0.05,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
} 
B <- seq(min(pairedend$score),max(pairedend$score),length.out = 100)
lines(getaxis(pairedend,B,xaxisfunction),getaxis(pairedend,B,yaxisfunction),type="l")
legend("bottomright",legend=S,col=m,lty=rep(1,10)) 

length(pairedend$pairedscore[pairedend$pairedscore > 0])

#Scan Bamsurgeon-trained model

colnames <- c("status","id","chrom","pos","score","depth","fracphased","fracstrand","frach","MAF","MAF_phased",
              "nC","fracC","Cavgbq","Cfrach","Cavgmapq","Cfracstrand","CMAF","CavgXL","CavgXV","CavgXR","CavgXB","Cavgpos","Cavgclip","Cavgind","Cfracclip","Cfracind","CavgASXS","Cavgpair",
              "Navgbq","Nfrach","Navgmapq","Nfracstrand","NMAF","NavgXL","NavgXV","NavgXR","NavgXB","Navgpos","Navgclip","Nfracclip","Nfracind","Navgind","NavgASXS","Navgpair",
              "weightedCbq","Cavgibq",
              "Mavgbq","Mavgmapq","Mfracstrand","Mavgpos","Mavgclip","Mavgind","Mfracclip","Mfracind","Mavgibq","weightedMbq","MavgASXS","Mavgpair",
              "Javgbq","Javgmapq","Jfracstrand","Javgpos","Javgclip","Javgind","Jfracclip","Jfracind","JavgASXS","Javgpair","MAFnorm","fracOtherAllele")

giabB <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/bamsurgeontrain/scan_norepeat.features_labeled.tsv",sep="\t",header=F,col.names=colnames)

plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="bamsurgeon training data - repeat filter only")
c <- 1
for (i in S) {
  t <- giabB[giabB$MAF >= i & giabB$MAF < i+0.05 & giabB$score > 0 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
} 
B <- seq(min(giabB$score),max(giabB$score),length.out = 100)
lines(getaxis(giabB[giabB$score > 0,],B,xaxisfunction),getaxis(giabB[giabB$score > 0,],B,yaxisfunction),type="l")
legend("bottomright",legend=S,col=m,lty=rep(1,10)) 

calls <- giabB
filteredB <- calls[calls$fracC>0 & calls$fracphased>0 & calls$Cfracclip<1 & calls$Mfracclip<1 & calls$Nfracclip<1 
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


plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="bamsurgeon training data -feature filtered")
c <- 1
for (i in S[2:10]) {
  t <- filteredB[filteredB$MAF >= i & filteredB$MAF < i+0.05 & filteredB$score > 0 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
} 
B <- seq(min(filteredB$score),max(filteredB$score),length.out = 100)
lines(getaxis(filteredB[filteredB$score > 0,],B,xaxisfunction),getaxis(filteredB[filteredB$score > 0,],B,yaxisfunction),type="l")
legend("bottomright",legend=S,col=m,lty=rep(1,10)) 

