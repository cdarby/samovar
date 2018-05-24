#MAF bins
S <- seq(0,0.45,by=0.05)

#### HCCT/N 1954 ####
#### In this section, method uses only its own region filter(s) ####
#### HapMuc ####

HapMuc <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/hapmuc_selffilter_labeled.tsv",sep="\t",header=F)
#V2 status; V22 Fisher; V25 hapmuc score; V26 description

#Calculate MAF from counts of number of reads
HapMuc$MAFT <- HapMuc$V9/(HapMuc$V8+HapMuc$V9)
#HapMuc$MAFN <- HapMuc$V11/(HapMuc$V10+HapMuc$V11)

#Fisher score
HapMuc$status <- HapMuc$V2[order(HapMuc$V22,decreasing = T)]
HapMuc$score <- HapMuc$V22[order(HapMuc$V22,decreasing = T)]

rank <- seq(1,length(HapMuc$V1))
HapMucFisherSelffilter <- cbind(sum(HapMuc$status[1:10]),sum(HapMuc$status[1:100]),sum(HapMuc$status[1:1000]),
                      sapply(S,function(i) length(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05])), #Number of calls in each bin
                      sapply(S,function(i) sum(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05])), #Number correct in each bin
                      sapply(S,function(i) sum(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                      sapply(S,function(i) length(HapMuc$MAFT[1:10][HapMuc$MAFT[1:10] >= i & HapMuc$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                      sapply(S,function(i) mean(HapMuc$score[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                      sapply(S,function(i) mean(rank[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                      sapply(S,function(i) sum(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                      sapply(S,function(i) length(HapMuc$MAFT[1:100][HapMuc$MAFT[1:100] >= i & HapMuc$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                      sapply(S,function(i) mean(HapMuc$score[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                      sapply(S,function(i) mean(rank[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                      sapply(S,function(i) sum(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                      sapply(S,function(i) length(HapMuc$MAFT[1:1000][HapMuc$MAFT[1:1000] >= i & HapMuc$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                      sapply(S,function(i) mean(HapMuc$score[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                      sapply(S,function(i) mean(rank[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(HapMucFisherSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/hapmuc_fisher_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F) 

#HapMuc score
H <- HapMuc[HapMuc$V25 != "-",]
s <- as.numeric(levels(H$V25))[H$V25]
H$score <- s[order(s,decreasing=T)]


rank <- seq(1,length(H$V1))
HapMucScoreSelffilter <- cbind(sum(H$status[1:10]),sum(H$status[1:100]),sum(H$status[1:1000]),
                    sapply(S,function(i) length(H$status[H$MAFT >= i & H$MAFT < i+0.05])), #Number of calls in each bin
                    sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05])), #Number correct in each bin
                    sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                    sapply(S,function(i) length(H$MAFT[1:10][H$MAFT[1:10] >= i & H$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                    sapply(S,function(i) mean(H$score[H$MAFT >= i & H$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                    sapply(S,function(i) mean(rank[H$MAFT >= i & H$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                    sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                    sapply(S,function(i) length(H$MAFT[1:100][H$MAFT[1:100] >= i & H$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                    sapply(S,function(i) mean(H$score[H$MAFT >= i & H$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                    sapply(S,function(i) mean(rank[H$MAFT >= i & H$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                    sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                    sapply(S,function(i) length(H$MAFT[1:1000][H$MAFT[1:1000] >= i & H$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                    sapply(S,function(i) mean(H$score[H$MAFT >= i & H$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                    sapply(S,function(i) mean(rank[H$MAFT >= i & H$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(HapMucScoreSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/hapmuc_score_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)






#### MosaicHunter Single ####

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/HCCT.mosaichunter.bed",sep="\t",header=F)
MHunter$MAFT <- MHunter$V15 / MHunter$V9
MHunter$status <- MHunter$V4
MHunter$score <- MHunter$V29
MHunter <- MHunter[order(MHunter$V29,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
MHSingleSelffilter <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
                     sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
                     sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
                     sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                     sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                     sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                     sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                     sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                     sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                     sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                     sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                     sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                     sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                     sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                     sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(MHSingleSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/mosaichunter_single_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)


#### MosaicHunter Pair ####

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/HCCPair.mosaichunter.bed",sep="\t",header=F)
MHunter$MAFT <- MHunter$V15 / MHunter$V9
MHunter$status <- MHunter$V4
MHunter$score <- MHunter$V50
MHunter <- MHunter[order(MHunter$V50,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
MHPairSelffilter <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(MH,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/mosaichunter_pair_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)



#### Samovar ####

Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/SamovarEdited.bed",sep="\t",header=F)
Samovar$MAFT <- Samovar$V18
Samovar$status <- Samovar$V4 / 2 #"real" is labeled 2; "bamsurgeon" is labeled 1
Samovar$score <- Samovar$V13
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
SmvSelffilter <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
            sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
            sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
            sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
            sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
            sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
            sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
            sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
            sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
            sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
            sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
            sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
            sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
            sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
            sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(SmvSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/samovar_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)


#### In this section, methods use the intersection of the region filters ####
#### HapMuc ####

HapMuc <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/HCCT.hapmuc.allvisible.tsv",sep="\t",header=F)
#V1 status; V22 Fisher; V25 hapmuc score; V26 description

#Calculate MAF from counts of number of reads
HapMuc$MAFT <- HapMuc$V9/(HapMuc$V8+HapMuc$V9)
HapMuc$MAFN <- HapMuc$V11/(HapMuc$V10+HapMuc$V11)

#Fisher score
HapMuc$status <- HapMuc$V1[order(HapMuc$V22,decreasing = T)]
HapMuc$score <- HapMuc$V22[order(HapMuc$V22,decreasing = T)]

#MAF bins
S <- seq(0,0.45,by=0.05)

rank <- seq(1,length(HapMuc$V1))
HapMucFisherAllvisible <- cbind(sum(HapMuc$status[1:10]),sum(HapMuc$status[1:100]),sum(HapMuc$status[1:1000]),
                      sapply(S,function(i) length(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05])), #Number of calls in each bin
                      sapply(S,function(i) sum(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05])), #Number correct in each bin
                      sapply(S,function(i) sum(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                      sapply(S,function(i) length(HapMuc$MAFT[1:10][HapMuc$MAFT[1:10] >= i & HapMuc$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                      sapply(S,function(i) mean(HapMuc$score[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                      sapply(S,function(i) mean(rank[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                      sapply(S,function(i) sum(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                      sapply(S,function(i) length(HapMuc$MAFT[1:100][HapMuc$MAFT[1:100] >= i & HapMuc$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                      sapply(S,function(i) mean(HapMuc$score[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                      sapply(S,function(i) mean(rank[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                      sapply(S,function(i) sum(HapMuc$status[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                      sapply(S,function(i) length(HapMuc$MAFT[1:1000][HapMuc$MAFT[1:1000] >= i & HapMuc$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                      sapply(S,function(i) mean(HapMuc$score[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                      sapply(S,function(i) mean(rank[HapMuc$MAFT >= i & HapMuc$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(HapMucFisher,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/hapmuc_fisher.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F) 

#HapMuc score
H <- HapMuc[HapMuc$V25 != "-",]
s <- as.numeric(levels(H$V25))[H$V25]
H$score <- s[order(s,decreasing=T)]


rank <- seq(1,length(H$V1))
HapMucScoreAllvisible <- cbind(sum(H$status[1:10]),sum(H$status[1:100]),sum(H$status[1:1000]),
                     sapply(S,function(i) length(H$status[H$MAFT >= i & H$MAFT < i+0.05])), #Number of calls in each bin
                     sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05])), #Number correct in each bin
                     sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                     sapply(S,function(i) length(H$MAFT[1:10][H$MAFT[1:10] >= i & H$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                     sapply(S,function(i) mean(H$score[H$MAFT >= i & H$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                     sapply(S,function(i) mean(rank[H$MAFT >= i & H$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                     sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                     sapply(S,function(i) length(H$MAFT[1:100][H$MAFT[1:100] >= i & H$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                     sapply(S,function(i) mean(H$score[H$MAFT >= i & H$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                     sapply(S,function(i) mean(rank[H$MAFT >= i & H$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                     sapply(S,function(i) sum(H$status[H$MAFT >= i & H$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                     sapply(S,function(i) length(H$MAFT[1:1000][H$MAFT[1:1000] >= i & H$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                     sapply(S,function(i) mean(H$score[H$MAFT >= i & H$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                     sapply(S,function(i) mean(rank[H$MAFT >= i & H$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(HapMucScoreAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/hapmuc_score.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F) 

#### MosaicHunter Single ####

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/HCCT.mosaichunter.allvisible.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V15 / MHunter$V9
MHunter$status <- MHunter$V4
MHunter$score <- MHunter$V29
MHunter <- MHunter[order(MHunter$V29,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
MHSingleAllvisible <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(MHSingleAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/mosaichunter_single.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)


#### MosaicHunter Pair ####

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/HCCPair.mosaichunter.allvisible.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V15 / MHunter$V9
MHunter$status <- MHunter$V4
MHunter$score <- MHunter$V50
MHunter <- MHunter[order(MHunter$V50,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
MHPairAllvisible <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(MHPairAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/mosaichunter_pair.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)




#### Samovar ####

Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/HCCT.editedsamovar.allvisible.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V18
Samovar$status <- Samovar$V4 / 2 #"real" is labeled 2; "bamsurgeon" is labeled 1
Samovar$score <- Samovar$V13
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
SmvAllvisible <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(SmvAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/samovar.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)





#### GIAB NA24385 ####
#### In this section, method uses only its own region filter(s) ####
#### MosaicHunter Single ####

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/MHunter/finalSingle.passed_labeled.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V12 / MHunter$V6
MHunter$status <- MHunter$V1
MHunter$score <- MHunter$V26 #Mosaic Posterior Probability
MHunter <- MHunter[order(MHunter$score,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
G_MHSingleSelffilter <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_MHSingleSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/mosaichunter_single_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/MHunter/finalSingleSub.passed_labeled.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V12 / MHunter$V6
MHunter$status <- MHunter$V1
MHunter$score <- MHunter$V26 #Mosaic Posterior Probability
MHunter <- MHunter[order(MHunter$score,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
G_MHSingleSubSelffilter <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
                              sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
                              sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
                              sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                              sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                              sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                              sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                              sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                              sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                              sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                              sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                              sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                              sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                              sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                              sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
write.table(G_MHSingleSubSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/mosaichunter_singleSub_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#### MosaicHunter Trio ####

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/MHunter/finalTrio.passed_labeled.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V12 / MHunter$V6
MHunter$status <- MHunter$V1
MHunter$score <- MHunter$V32 #Mosaic Posterior Probability (child)
MHunter <- MHunter[order(MHunter$score,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
G_MHTrioSelffilter <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_MHTrioSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/mosaichunter_pair_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/MHunter/finalTrioSub.passed_labeled.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V12 / MHunter$V6
MHunter$status <- MHunter$V1
MHunter$score <- MHunter$V32 #Mosaic Posterior Probability (child)
MHunter <- MHunter[order(MHunter$score,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
G_MHTrioSubSelffilter <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
                            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
                            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
                            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_MHTrioSubSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/mosaichunter_pairSub_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#### Samovar ####

#Samovar
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/editedsitestrain/filteredSamovar_linkage.labeled.tsv",sep="\t",header=F)
Samovar <- Samovar[Samovar$V5 >= 0.05,] #linkage filter
Samovar$MAFT <- Samovar$V15
Samovar$status <- Samovar$V1
Samovar$score <- Samovar$V13 #samovar score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvSelffilter <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#Paired
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/editedsitestrain/filteredPair.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V9
Samovar$status <- Samovar$V1
Samovar$score <- Samovar$V6 #paired score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvPairSelffilter <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvPairSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_paired_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#No-Phasing
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/editedsitestrain/filteredUnphased.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V9
Samovar$status <- Samovar$V1
Samovar$score <- Samovar$V5 #unphased score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvUnphasedSelffilter <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvUnphasedSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_unphased_selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)


#Samovar - Subsample
#Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/NEWfilteredSamovarSubsample_linkage.labeled.tsv",sep="\t",header=F)
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/editedsitestrain/subsample/filteredSamovarSubsample_linkage.labeled.tsv",sep="\t",header=F)
Samovar <- Samovar[Samovar$V5 >= 0.05,] #linkage filter
Samovar$MAFT <- Samovar$V15
Samovar$status <- Samovar$V1
Samovar$score <- Samovar$V13 #samovar score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvSubSelffilter <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvSubSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/NEWsamovar_subsample.selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#Paired - Subsample
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/editedsitestrain/subsample/filteredPairSubsample.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V12
Samovar$status <- Samovar$V1
Samovar$score <- Samovar$V6 #paired score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvSubPairSelffilter <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvSubPairSelffilter,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_subsample_paired.selffilter.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#No-Phasing - Subsample
#TODO check if correct file
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giab/editedsitestrain/subsample/filteredUnphasedSubsample.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V12
Samovar$status <- Samovar$V1
Samovar$score <- Samovar$V5 #unphased score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvSubUnphasedSelffilter <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(Smv,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_subsample_unphased.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)




#### In this section, methods use the intersection of the region filters ####
#### MosaicHunter Single ####

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/finalSingle.passed_labeled1.tsv.allvisible.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V15 / MHunter$V9
MHunter$status <- MHunter$V4
MHunter$score <- MHunter$V29 #Mosaic Posterior Probability
MHunter <- MHunter[order(MHunter$score,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
G_MHSingleAllvisible <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_MHSingleAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/mosaichunter_single_allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/finalSingleSub.passed_labeled.bed.allvisible.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V15 / MHunter$V9
MHunter$status <- MHunter$V4
MHunter$score <- MHunter$V29 #Mosaic Posterior Probability
MHunter <- MHunter[order(MHunter$score,decreasing=T),]

rank <- seq(1,length(MHunter$V1))
G_MHSingleSubAllvisible <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
                              sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
                              sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
                              sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                              sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                              sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                              sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                              sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                              sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                              sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                              sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                              sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                              sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                              sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                              sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_MHSingleSubAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/mosaichunterSub_single_allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)


#### MosaicHunter Trio ####

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/finalTrio.passed_labeled1.tsv.allvisible.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V15 / MHunter$V9
MHunter$status <- MHunter$V4
MHunter$score <- MHunter$V35 #Mosaic Posterior Probability (child)
MHunter <- MHunter[order(MHunter$score,decreasing=T),]


rank <- seq(1,length(MHunter$V1))
G_MHTrioAllvisible <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_MHTrioAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/mosaichunter_trio_allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

MHunter <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/finalTrioSub.passed_labeled.bed.allvisible.tsv",sep="\t",header=F)
MHunter$MAFT <- MHunter$V15 / MHunter$V9
MHunter$status <- MHunter$V4
MHunter$score <- MHunter$V35 #Mosaic Posterior Probability (child)
MHunter <- MHunter[order(MHunter$score,decreasing=T),]


rank <- seq(1,length(MHunter$V1))
G_MHTrioSubAllvisible <- cbind(sum(MHunter$status[1:10]),sum(MHunter$status[1:100]),sum(MHunter$status[1:1000]),
                            sapply(S,function(i) length(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number of calls in each bin
                            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05])), #Number correct in each bin
                            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
                            sapply(S,function(i) length(MHunter$MAFT[1:10][MHunter$MAFT[1:10] >= i & MHunter$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
                            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
                            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
                            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
                            sapply(S,function(i) length(MHunter$MAFT[1:100][MHunter$MAFT[1:100] >= i & MHunter$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
                            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
                            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
                            sapply(S,function(i) sum(MHunter$status[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
                            sapply(S,function(i) length(MHunter$MAFT[1:1000][MHunter$MAFT[1:1000] >= i & MHunter$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
                            sapply(S,function(i) mean(MHunter$score[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
                            sapply(S,function(i) mean(rank[MHunter$MAFT >= i & MHunter$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
write.table(G_MHTrioSubAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/HCCTAnalysis/mosaichunter_trioSub_allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)


#### Samovar ####

#Samovar
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/filteredSamovar_linkage1.labeled.tsv.allvisible.tsv",sep="\t",header=F)
Samovar <- Samovar[Samovar$V8 >= 0.05,] #linkage filter
Samovar$MAFT <- Samovar$V18
Samovar$status <- Samovar$V4
Samovar$score <- Samovar$V16 #samovar score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvAllvisible <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#Paired
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/filteredPair1.tsv.allvisible.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V12
Samovar$status <- Samovar$V4
Samovar$score <- Samovar$V9 #paired score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvPairAllvisible <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvPairAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_paired.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#No-Phasing
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/filteredUnphased1.tsv.allvisible.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V12
Samovar$status <- Samovar$V4
Samovar$score <- Samovar$V8 #unphased score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvUnphasedAllvisible <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvUnphasedAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_unphased.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)


#Samovar - Subsample
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/filteredSamovarSubsample_linkage1.labeled.tsv.allvisible.tsv",sep="\t",header=F)
Samovar <- Samovar[Samovar$V8 >= 0.05,] #linkage filter
Samovar$MAFT <- Samovar$V18
Samovar$status <- Samovar$V4
Samovar$score <- Samovar$V16 #samovar score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvSubAllvisible <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvSubAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_subsample.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#Paired - Subsample
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/filteredPairSubsample1.tsv.allvisible.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V12
Samovar$status <- Samovar$V4
Samovar$score <- Samovar$V9 #paired score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvSubPairAllvisible <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvSubPairAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_subsample_paired.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)

#No-Phasing - Subsample
Samovar <- read.delim("~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/filteredUnphasedSubsample1.tsv.allvisible.tsv",sep="\t",header=F)
Samovar$MAFT <- Samovar$V12
Samovar$status <- Samovar$V4
Samovar$score <- Samovar$V8 #unphased score
Samovar <- Samovar[order(Samovar$score,decreasing=T),]

rank <- seq(1,length(Samovar$V1))
G_SmvSubUnphasedAllvisible <- cbind(sum(Samovar$status[1:10]),sum(Samovar$status[1:100]),sum(Samovar$status[1:1000]),
             sapply(S,function(i) length(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number of calls in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05])), #Number correct in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Number correct of top 10 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:10][Samovar$MAFT[1:10] >= i & Samovar$MAFT[1:10] < i+0.05])), #Number in overall top 10 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:10])), #Mean score of top 10 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean rank of top 10 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Number correct of top 100 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:100][Samovar$MAFT[1:100] >= i & Samovar$MAFT[1:100] < i+0.05])), #Number in overall top 100 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:100])), #Mean score of top 100 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean rank of top 100 in each bin
             sapply(S,function(i) sum(Samovar$status[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Number correct of top 1000 in each bin
             sapply(S,function(i) length(Samovar$MAFT[1:1000][Samovar$MAFT[1:1000] >= i & Samovar$MAFT[1:1000] < i+0.05])), #Number in overall top 1000 from each bin
             sapply(S,function(i) mean(Samovar$score[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000])), #Mean score of top 1000 in each bin
             sapply(S,function(i) mean(rank[Samovar$MAFT >= i & Samovar$MAFT < i+0.05][1:1000]))) #Mean rank of top 1000 in each bin
#write.table(G_SmvSubUnphasedAllvisible,"~/Documents/lab/samovar/testdata_DO_NOT_COMMIT/giabAnalysis/samovar_subsample_unphased.allvisible.table.tsv",row.names=F,col.names=F,sep="\t",quote=F)




#### Plots HCC1954 ####
par(mfrow=c(2,3))

#Top 10 - count correct
plot(seq(0,0.45,0.05),HapMucFisherSelffilter[1:10,6],ylim=c(0,10),type="b",pch=19,col="green",main="Count correct of top 10 (HCC1954)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0,0.45,0.05),HapMucScoreSelffilter[1:10,6],type="b",pch=19,col="forestgreen")
points(seq(0,0.45,0.05),MHPairSelffilter[1:10,6],type="b",pch=19,col="cyan")
points(seq(0,0.45,0.05),MHSingleSelffilter[1:10,6],type="b",pch=19,col="blue")
points(seq(0,0.45,0.05),SmvSelffilter[1:10,6],type="b",pch=19,col="red")

axis(1,at=seq(0,0.45,0.05))
axis(2,at=seq(0,10))
legend("bottomright",legend=c(paste("HapMucFisher",HapMucFisherSelffilter[1]),paste("HapMucScore",HapMucScoreSelffilter[1]),
                              paste("MosaicHunterPair",MHPairSelffilter[1]),paste("MosaicHunterSingle",MHSingleSelffilter[1]),
                              paste("Samovar",SmvSelffilter[1])),col=c("green","forestgreen","cyan","blue","red"),pch=19)

#Top 100 - count correct
plot(seq(0,0.45,0.05),HapMucFisherSelffilter[1:10,10],ylim=c(0,100),type="b",pch=19,col="green",main="Count correct of top 100 (HCC1954)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0,0.45,0.05),HapMucScoreSelffilter[1:10,10],type="b",pch=19,col="forestgreen")
points(seq(0,0.45,0.05),MHPairSelffilter[1:10,10],type="b",pch=19,col="cyan")
points(seq(0,0.45,0.05),MHSingleSelffilter[1:10,10],type="b",pch=19,col="blue")
points(seq(0,0.45,0.05),SmvSelffilter[1:10,10],type="b",pch=19,col="red")

axis(1,at=seq(0,0.45,0.05))
axis(2,at=seq(0,100,5))
legend("bottomright",legend=c(paste("HapMucFisher",HapMucFisherSelffilter[1,2]),paste("HapMucScore",HapMucScoreSelffilter[1,2]),
                              paste("MosaicHunterPair",MHPairSelffilter[1,2]),paste("MosaicHunterSingle",MHSingleSelffilter[1,2]),
                              paste("Samovar",SmvSelffilter[1,2])),col=c("green","forestgreen","cyan","blue","red"),pch=19)

#Accuracy of all calls
plot(seq(0,0.45,0.05),HapMucFisherSelffilter[1:10,5]/HapMucFisherSelffilter[1:10,4],ylim=c(0,1),type="b",pch=19,col="green",main="Fraction correct of all calls made (HCC1954)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0,0.45,0.05),HapMucScoreSelffilter[1:10,5]/HapMucScoreSelffilter[1:10,4],type="b",pch=19,col="forestgreen")
points(seq(0,0.45,0.05),MHPairSelffilter[1:10,5]/MHPairSelffilter[1:10,4],type="b",pch=19,col="cyan")
points(seq(0,0.45,0.05),MHSingleSelffilter[1:10,5]/MHSingleSelffilter[1:10,4],type="b",pch=19,col="blue")
points(seq(0,0.45,0.05),SmvSelffilter[1:10,5]/SmvSelffilter[1:10,4],type="b",pch=19,col="red")

axis(1,at=seq(0,0.45,0.05))
axis(2,at=seq(0,1,0.1))
legend("topleft",legend=c(paste("HapMucFisher",sum(HapMucFisherSelffilter[1:10,4])),paste("HapMucScore",sum(HapMucScoreSelffilter[1:10,4])),
                              paste("MosaicHunterPair",sum(MHPairSelffilter[1:10,4])),paste("MosaicHunterSingle",sum(MHSingleSelffilter[1:10,4])),
                              paste("Samovar",sum(SmvSelffilter[1:10,4]))),col=c("green","forestgreen","cyan","blue","red"),pch=19)




#Top 10 - count correct
plot(seq(0,0.45,0.05),HapMucFisherAllvisible[1:10,6],ylim=c(0,10),type="b",pch=19,col="green",main="Count correct of top 10 (HCC1954)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0,0.45,0.05),HapMucScoreAllvisible[1:10,6],type="b",pch=19,col="forestgreen")
points(seq(0,0.45,0.05),MHPairAllvisible[1:10,6],type="b",pch=19,col="cyan")
points(seq(0,0.45,0.05),MHSingleAllvisible[1:10,6],type="b",pch=19,col="blue")
points(seq(0,0.45,0.05),SmvAllvisible[1:10,6],type="b",pch=19,col="red")

axis(1,at=seq(0,0.45,0.05))
axis(2,at=seq(0,10))
legend("bottomright",legend=c(paste("HapMucFisher",HapMucFisherAllvisible[1]),paste("HapMucScore",HapMucScoreAllvisible[1]),
                              paste("MosaicHunterPair",MHPairAllvisible[1]),paste("MosaicHunterSingle",MHSingleAllvisible[1]),
                              paste("Samovar",SmvAllvisible[1])),col=c("green","forestgreen","cyan","blue","red"),pch=19)

#Top 100 - count correct
plot(seq(0,0.45,0.05),HapMucFisherAllvisible[1:10,10],ylim=c(0,100),type="b",pch=19,col="green",main="Count correct of top 100 (HCC1954)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0,0.45,0.05),HapMucScoreAllvisible[1:10,10],type="b",pch=19,col="forestgreen")
points(seq(0,0.45,0.05),MHPairAllvisible[1:10,10],type="b",pch=19,col="cyan")
points(seq(0,0.45,0.05),MHSingleAllvisible[1:10,10],type="b",pch=19,col="blue")
points(seq(0,0.45,0.05),SmvAllvisible[1:10,10],type="b",pch=19,col="red")

axis(1,at=seq(0,0.45,0.05))
axis(2,at=seq(0,100,5))
legend("bottomright",legend=c(paste("HapMucFisher",HapMucFisherAllvisible[1,2]),paste("HapMucScore",HapMucScoreAllvisible[1,2]),
                              paste("MosaicHunterPair",MHPairAllvisible[1,2]),paste("MosaicHunterSingle",MHSingleAllvisible[1,2]),
                              paste("Samovar",SmvAllvisible[1,2])),col=c("green","forestgreen","cyan","blue","red"),pch=19)

#Accuracy of all calls
plot(seq(0,0.45,0.05),HapMucFisherAllvisible[1:10,5]/HapMucFisherAllvisible[1:10,4],ylim=c(0,1),type="b",pch=19,col="green",main="Fraction correct of all calls made (HCC1954)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0,0.45,0.05),HapMucScoreAllvisible[1:10,5]/HapMucScoreAllvisible[1:10,4],type="b",pch=19,col="forestgreen")
points(seq(0,0.45,0.05),MHPairAllvisible[1:10,5]/MHPairAllvisible[1:10,4],type="b",pch=19,col="cyan")
points(seq(0,0.45,0.05),MHSingleAllvisible[1:10,5]/MHSingleAllvisible[1:10,4],type="b",pch=19,col="blue")
points(seq(0,0.45,0.05),SmvAllvisible[1:10,5]/SmvAllvisible[1:10,4],type="b",pch=19,col="red")

axis(1,at=seq(0,0.45,0.05))
axis(2,at=seq(0,1,0.1))
legend("topleft",legend=c(paste("HapMucFisher",sum(HapMucFisherAllvisible[1:10,4])),paste("HapMucScore",sum(HapMucScoreAllvisible[1:10,4])),
                          paste("MosaicHunterPair",sum(MHPairAllvisible[1:10,4])),paste("MosaicHunterSingle",sum(MHSingleAllvisible[1:10,4])),
                          paste("Samovar",sum(SmvAllvisible[1:10,4]))),col=c("green","forestgreen","cyan","blue","red"),pch=19)


#### Plots GIAB NA24385 ####

#Top 10 - count correct
plot(seq(0.05,0.45,0.05),G_SmvPairSelffilter[2:10,6],ylim=c(0,10),type="b",pch=19,col="maroon",main="Count correct of top 10 (bamsurgeon of NA24385)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvUnphasedSelffilter[2:10,6],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSelffilter[2:10,6],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSelffilter[2:10,6],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSelffilter[2:10,6],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,10))
legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvPairSelffilter[1]),paste("SamovarUnphased",G_SmvUnphasedSelffilter[1]),
                              paste("Samovar",G_SmvSelffilter[1]), paste("MosaicHunterTrio",G_MHTrioSelffilter[1]),paste("MosaicHunterSingle",G_MHSingleSelffilter[1])),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)

#Top 100 - count correct
plot(seq(0.05,0.45,0.05),G_SmvPairSelffilter[2:10,10],ylim=c(0,100),type="b",pch=19,col="maroon",main="Count correct of top 100 (bamsurgeon of NA24385)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvUnphasedSelffilter[2:10,10],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSelffilter[2:10,10],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSelffilter[2:10,10],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSelffilter[2:10,10],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,100,10))
legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvPairSelffilter[1,2]),paste("SamovarUnphased",G_SmvUnphasedSelffilter[1,2]),
                              paste("Samovar",G_SmvSelffilter[1,2]), paste("MosaicHunterTrio",G_MHTrioSelffilter[1,2]),paste("MosaicHunterSingle",G_MHSingleSelffilter[1,2])),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)

#Accuracy of all calls
plot(seq(0.05,0.45,0.05),G_SmvPairSelffilter[2:10,5]/G_SmvPairSelffilter[2:10,4],ylim=c(0,1),type="b",pch=19,col="maroon",main="Fraction correct of all calls made (bamsurgeon of NA24385)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvUnphasedSelffilter[2:10,5]/G_SmvUnphasedSelffilter[2:10,4],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSelffilter[2:10,5]/G_SmvSelffilter[2:10,4],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSelffilter[2:10,5]/G_MHTrioSelffilter[2:10,4],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSelffilter[2:10,5]/G_MHSingleSelffilter[2:10,4],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,1,0.1))
legend("bottomright",legend=c(paste("SamovarShortRead",sum(G_SmvPairSelffilter[1:10,4])),paste("SamovarUnphased",sum(G_SmvUnphasedSelffilter[1:10,4])),
                              paste("Samovar",sum(G_SmvSelffilter[1:10,4])), paste("MosaicHunterTrio",sum(G_MHTrioSelffilter[1:10,4])),paste("MosaicHunterSingle",sum(G_MHSingleSelffilter[1:10,4]))),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)



#Top 10 - count correct
plot(seq(0.05,0.45,0.05),G_SmvPairAllvisible[2:10,6],ylim=c(0,10),type="b",pch=19,col="maroon",main="Count correct of top 10 (bamsurgeon of NA24385)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvUnphasedAllvisible[2:10,6],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvAllvisible[2:10,6],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioAllvisible[2:10,6],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleAllvisible[2:10,6],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,10))
legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvPairAllvisible[1]),paste("SamovarUnphased",G_SmvUnphasedAllvisible[1]),
                              paste("Samovar",G_SmvAllvisible[1]), paste("MosaicHunterTrio",G_MHTrioAllvisible[1]),paste("MosaicHunterSingle",G_MHSingleAllvisible[1])),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)

#Top 100 - count correct
plot(seq(0.05,0.45,0.05),G_SmvPairAllvisible[2:10,10],ylim=c(0,100),type="b",pch=19,col="maroon",main="Count correct of top 100 (bamsurgeon of NA24385)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvUnphasedAllvisible[2:10,10],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvAllvisible[2:10,10],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioAllvisible[2:10,10],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleAllvisible[2:10,10],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,100,10))
legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvPairAllvisible[1,2]),paste("SamovarUnphased",G_SmvUnphasedAllvisible[1,2]),
                              paste("Samovar",G_SmvAllvisible[1,2]), paste("MosaicHunterTrio",G_MHTrioAllvisible[1,2]),paste("MosaicHunterSingle",G_MHSingleAllvisible[1,2])),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)

#Accuracy of all calls
plot(seq(0.05,0.45,0.05),G_SmvPairAllvisible[2:10,5]/G_SmvPairAllvisible[2:10,4],ylim=c(0,1),type="b",pch=19,col="maroon",main="Fraction correct of all calls made (bamsurgeon of NA24385)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvUnphasedAllvisible[2:10,5]/G_SmvUnphasedAllvisible[2:10,4],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvAllvisible[2:10,5]/G_SmvAllvisible[2:10,4],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioAllvisible[2:10,5]/G_MHTrioAllvisible[2:10,4],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleAllvisible[2:10,5]/G_MHSingleAllvisible[2:10,4],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,1,0.1))
legend("bottomright",legend=c(paste("SamovarShortRead",sum(G_SmvPairAllvisible[1:10,4])),paste("SamovarUnphased",sum(G_SmvUnphasedAllvisible[1:10,4])),
                              paste("Samovar",sum(G_SmvAllvisible[1:10,4])), paste("MosaicHunterTrio",sum(G_MHTrioAllvisible[1:10,4])),paste("MosaicHunterSingle",sum(G_MHSingleAllvisible[1:10,4]))),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)



#Top 10 - count correct
plot(seq(0.05,0.45,0.05),G_SmvSubPairSelffilter[2:10,6],ylim=c(0,10),type="b",pch=19,col="maroon",main="Count correct of top 10 (subsampled; bamsurgeon of NA24385)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvSubUnphasedSelffilter[2:10,6],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSubSelffilter[2:10,6],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSubSelffilter[2:10,6],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSubSelffilter[2:10,6],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,10))
legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvSubPairSelffilter[1]),paste("SamovarUnphased",G_SmvSubUnphasedSelffilter[1]),
                              paste("Samovar",G_SmvSubSelffilter[1]), paste("MosaicHunterTrio",G_MHTrioSubSelffilter[1]),paste("MosaicHunterSingle",G_MHSingleSubSelffilter[1])),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)

#Top 100 - count correct
plot(seq(0.05,0.45,0.05),G_SmvSubPairSelffilter[2:10,10],ylim=c(0,100),type="b",pch=19,col="maroon",main="Count correct of top 100 (subsampled; bamsurgeon of NA24385)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvSubUnphasedSelffilter[2:10,10],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSubSelffilter[2:10,10],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSubSelffilter[2:10,10],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSubSelffilter[2:10,10],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,100,10))
legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvSubPairSelffilter[1,2]),paste("SamovarUnphased",G_SmvSubUnphasedSelffilter[1,2]),
                              paste("Samovar",G_SmvSubSelffilter[1,2]), paste("MosaicHunterTrio",G_MHTrioSubSelffilter[1,2]),paste("MosaicHunterSingle",G_MHSingleSubSelffilter[1,2])),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)

#Accuracy of all calls
plot(seq(0.05,0.45,0.05),G_SmvSubPairSelffilter[2:10,5]/G_SmvSubPairSelffilter[2:10,4],ylim=c(0,1),type="b",pch=19,col="maroon",main="Fraction correct of all calls made (subsampled; bamsurgeon of NA24385)\nIndividual Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvSubUnphasedSelffilter[2:10,5]/G_SmvSubUnphasedSelffilter[2:10,4],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSubSelffilter[2:10,5]/G_SmvSubSelffilter[2:10,4],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSubSelffilter[2:10,5]/G_MHTrioSubSelffilter[2:10,4],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSubSelffilter[2:10,5]/G_MHSingleSubSelffilter[2:10,4],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,1,0.1))
legend("bottomright",legend=c(paste("SamovarShortRead",sum(G_SmvSubPairSelffilter[1:10,4])),paste("SamovarUnphased",sum(G_SmvSubUnphasedSelffilter[1:10,4])),
                              paste("Samovar",sum(G_SmvSubSelffilter[1:10,4])), paste("MosaicHunterTrio",sum(G_MHTrioSubSelffilter[1:10,4])),paste("MosaicHunterSingle",sum(G_MHSingleSubSelffilter[1:10,4]))),
       col=c("maroon","darkorange2","red","cyan","blue"),pch=19)



#Top 10 - count correct
plot(seq(0.05,0.45,0.05),G_SmvSubPairAllvisible[2:10,6],ylim=c(0,10),type="b",pch=19,col="maroon",main="Count correct of top 10 (subsampled; bamsurgeon of NA24385)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvSubUnphasedAllvisible[2:10,6],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSubAllvisible[2:10,6],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSubAllvisible[2:10,6],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSubAllvisible[2:10,6],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,10))
legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvSubPairAllvisible[1]),paste("SamovarUnphased",G_SmvSubUnphasedAllvisible[1]),
                              paste("Samovar",G_SmvSubAllvisible[1]), paste("MosaicHunterTrio",G_MHTrioSubAllvisible[1]),paste("MosaicHunterSingle",G_MHSingleSubAllvisible[1])),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)
# legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvSubPairAllvisible[1]),paste("SamovarUnphased",G_SmvSubUnphasedAllvisible[1]),
#                               paste("Samovar",G_SmvSubAllvisible[1]), paste("MosaicHunterTrio"),paste("MosaicHunterSingle")),
#        col=c("maroon","darkorange2","red","cyan","blue"),pch=19)

#Top 100 - count correct
plot(seq(0.05,0.45,0.05),G_SmvSubPairAllvisible[2:10,10],ylim=c(0,100),type="b",pch=19,col="maroon",main="Count correct of top 100 (subsampled; bamsurgeon of NA24385)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvSubUnphasedAllvisible[2:10,10],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSubAllvisible[2:10,10],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSubAllvisible[2:10,10],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSubAllvisible[2:10,10],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,100,10))
legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvSubPairAllvisible[1,2]),paste("SamovarUnphased",G_SmvSubUnphasedAllvisible[1,2]),
                              paste("Samovar",G_SmvSubAllvisible[1,2]), paste("MosaicHunterTrio",G_MHTrioSubAllvisible[1,2]),paste("MosaicHunterSingle",G_MHSingleSubAllvisible[1,2])),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)
# legend("bottomright",legend=c(paste("SamovarShortRead",G_SmvSubPairAllvisible[1,2]),paste("SamovarUnphased",G_SmvSubUnphasedAllvisible[1,2]),
#                               paste("Samovar",G_SmvSubAllvisible[1,2]), paste("MosaicHunterTrio"),paste("MosaicHunterSingle")),
#                               col=c("maroon","darkorange2","red","cyan","blue"),pch=19)

#Accuracy of all calls
plot(seq(0.05,0.45,0.05),G_SmvSubPairAllvisible[2:10,5]/G_SmvSubPairAllvisible[2:10,4],ylim=c(0,1),type="b",pch=19,col="maroon",main="Fraction correct of all calls made (subsampled; bamsurgeon of NA24385)\nCombined Filters",axes=F,xlab="MAF (Bin)",ylab="Number correct")
points(seq(0.05,0.45,0.05),G_SmvSubUnphasedAllvisible[2:10,5]/G_SmvSubUnphasedAllvisible[2:10,4],type="b",pch=19,col="darkorange2")
points(seq(0.05,0.45,0.05),G_SmvSubAllvisible[2:10,5]/G_SmvSubAllvisible[2:10,4],type="b",pch=19,col="red")
points(seq(0.05,0.45,0.05),G_MHTrioSubAllvisible[2:10,5]/G_MHTrioSubAllvisible[2:10,4],type="b",pch=19,col="cyan")
points(seq(0.05,0.45,0.05),G_MHSingleSubAllvisible[2:10,5]/G_MHSingleSubAllvisible[2:10,4],type="b",pch=19,col="blue")

axis(1,at=seq(0.05,0.45,0.05))
axis(2,at=seq(0,1,0.1))
legend("bottomright",legend=c(paste("SamovarShortRead",sum(G_SmvSubPairAllvisible[1:10,4])),paste("SamovarUnphased",sum(G_SmvSubUnphasedAllvisible[1:10,4])),
                              paste("Samovar",sum(G_SmvSubAllvisible[1:10,4])), paste("MosaicHunterTrio",sum(G_MHTrioSubAllvisible[1:10,4])),paste("MosaicHunterSingle",sum(G_MHSingleSubAllvisible[1:10,4]))),
                              col=c("maroon","darkorange2","red","cyan","blue"),pch=19)
# legend("bottomright",legend=c(paste("SamovarShortRead",sum(G_SmvSubPairAllvisible[1:10,4])),paste("SamovarUnphased",sum(G_SmvSubUnphasedAllvisible[1:10,4])),
#                               paste("Samovar",sum(G_SmvSubAllvisible[1:10,4])), paste("MosaicHunterTrio"),paste("MosaicHunterSingle")),
#                               col=c("maroon","darkorange2","red","cyan","blue"),pch=19)


