#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

library(flux)

#FPR,1-Specificity = false positive / condition negative
#TPR,Sensitivity,Recall = true positive / condition positive

#0 = germline; 1 = mosaic (simulated)

fprfromthreshold <- function(p,A) { #predicted positives; all
  return(length(p$status[p$status==0])/(length(A$status[A$status==01])))
}

tprfromthreshold <- function(p,A) { #predicted positives; all
  return(length(p$status[p$status==1])/(length(A$status[A$status==1])))
}

getaxisnumericalth <- function(sortedvector,steps,partitionfunction){
  theaxis <- sapply(steps,function(i) partitionfunction(sortedvector[sortedvector$score >= i,],sortedvector)) #numerical threshold
  min0 <- min(theaxis)
  normalize <- max(theaxis)-min(theaxis)
  return(sapply(theaxis,function(i) (i-min0)/normalize))
}

getaxis <- getaxisnumericalth
A<-seq(0,1,0.01)
xaxisfunction <- fprfromthreshold
yaxisfunction <- tprfromthreshold


#X <- read.delim(header=F,"/Users/charlotte/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/depth_mosaicsim/allfeatures/3000.txt",sep="\t",col.names=c("status","score"))
X <- read.delim(header=F,args[1],sep="\t",col.names=c("status","score"))
#X <- X[order(X$score,decreasing=T),]

#plot(getaxis(X,A,xaxisfunction),getaxis(X,A,yaxisfunction),type="l",xlab="FPR: false positive / condition negative",ylab="TPR: true positive / condition positive",xlim=c(0,1),ylim=c(0,1))
print(auc(getaxis(X,A,xaxisfunction),getaxis(X,A,yaxisfunction)))
