HapMuc <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCPair/hapmuc/calls2.txt",sep="\t",header=F)
#V1 status; V22 Fisher; V25 hapmuc score; V26 description

HapMuc$MAFT <- HapMuc$V9/(HapMuc$V8+HapMuc$V9)
HapMuc$MAFN <- HapMuc$V11/(HapMuc$V10+HapMuc$V11)

HapMuc$status <- HapMuc$V1[order(HapMuc$V22,decreasing = T)]
HapMuc$score <- HapMuc$V22[order(HapMuc$V22,decreasing = T)]
B <- seq(0,12.7,length.out = 100)
plot(getaxis(HapMuc,B,xaxisfunction),getaxis(HapMuc,B,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
auc(getaxis(HapMuc,B,xaxisfunction),getaxis(HapMuc,B,yaxisfunction))
#0.892 - By fisher score, which every site has (15484 true positive).

HapMuc <- HapMuc[order(HapMuc$score,decreasing=TRUE),]

par(mfrow=c(4,3))
for (i in seq(0,0.45,by=0.05)) {
  t <- HapMuc[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05 ,]
  B <- seq(min(t$score),max(t$score),length.out = 100)
  plot(main=paste(i,length(t$status),round(auc(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction)),5)),getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
}




H <- HapMuc[HapMuc$V25 != "-",]
s <- as.numeric(levels(H$V25))[H$V25]
H$score <- s[order(s,decreasing=T)]
#H$status <- H$status[order(H$score,decreasing=T)]

B <- seq(-86,67,length.out = 500)
B <- seq(1,67,length.out = 200)

plot(getaxis(H,B,xaxisfunction),getaxis(H,B,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
auc(getaxis(H,B,xaxisfunction),getaxis(H,B,yaxisfunction))
# 0.557 - By hapmuc score, which 36774 have (4263 true positive). (Don't consider sites with no hapmuc score)
# 0.928 - By hapmuc score if > 0 , which 7398 have (2075 true positive). (Don't consider sites with no hapmuc score)
# 0.924 - By hapmuc score if > 1 , which 6084 have (2043 true positive). (Don't consider sites with no hapmuc score)

par(mfrow=c(4,3))
for (i in seq(0,0.45,by=0.05)) {
  t <- H[H$MAFT >= i & H$MAFT < i+0.05 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  plot(main=paste(i,length(t$status), round(auc(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction)),5)),getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
}

K <- colorRampPalette(c("red","purple","cyan"))
m <- K(10)
S <- seq(0,0.45,by=0.05)

plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="HapMuc (hapmuc score ranked)")
c <- 1
for (i in S) {
  t <- H[H$MAFT >= i & H$MAFT < i+0.05 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
}
legend("bottomright",legend=S,col=m,lty=rep(1,10)) 


plot(c(0,1),c(0,1),type="l",lty=2,xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",main="HapMuc (fisher score ranked)")
c <- 1
for (i in S) {
  t <- HapMuc[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05 ,] #& H$score > 0
  B <- seq(min(t$score),max(t$score),length.out = 100)
  lines(getaxis(t,B,xaxisfunction),getaxis(t,B,yaxisfunction),col=m[c])
  c <- c+1
}
legend("bottomright",legend=S,col=m,lty=rep(1,10)) 


rank <- seq(1,159979)
rank <- seq(1,36774)

sapply(S,function(i) length(H$status[H$MAFT >= i & H$MAFT < i+0.05]))
sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05]))

sum(H$status[1:10])
sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05][1:10]))
sapply(S,function(i) length(H$MAFT[1:10][H$MAFT[1:10] >= i & H$MAFT[1:10] < i+0.05]))
sapply(S,function(i) mean(H$score[H$MAFT >= i & H$MAFT < i+0.05][1:10]))
sapply(S,function(i) mean(rank[H$MAFT >= i & H$MAFT < i+0.05][1:100]))








