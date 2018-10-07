setwd("D:/work/pkt206")
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)
library("BiocParallel")
library(edgeR)
# p53 wt untr vs p53 wt tr
wunt1<-read.table(file="countTableOfp53_wt_untr.txt",header=TRUE, sep="\t")
wt1<-read.table(file="countTableOfp53_wt_tr.txt",header=TRUE, sep="\t")
countdata<-cbind(wunt1,wt1)
group1<-factor(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))

library(edgeR)
y<-DGEList(counts=countdata,group=group1)
#perform the TMM normalization
y<-calcNormFactors(y)
design<-model.matrix(~group1)
#To estimate common dispersion and tagwise dispersions in one run
y<-estimateDisp(y,design)
#et<-exactTest(y)
#topTags(et)
#use GLM  likelihood ratio test to determing Differential expressed genes
fit<-glmFit(y,design)
# compare mu vs wt
lrt.2vs1 <- glmLRT(fit, coef=2)

result_step_1<-topTags(lrt.2vs1,n=nrow(lrt.2vs1),adjust.method="BH", sort.by="logFC")

# the output is a top table sorted by logFC of each transcript
#write.table(result_step_1,file="GLMResultsOf4SamplesbylogFC.txt",sep="\t",quote=FALSE)

#annotation
library("org.Hs.eg.db")
d0<-result_step_1$table
d1<-cbind(rownames(d0),d0)
colnames(d1)<-c("Id","logFC","unshrunk.logFC","logCPM","PValue","FDR")

d1$symbol <- mapIds(org.Hs.eg.db,
                    keys=as.character(d1$Id),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")

d2<-na.omit(d1)
write.csv(d2,file="p53_wt_tr_vs_untr.csv",sep="\t",quote=FALSE)


rd<-unique(as.character(d2$symbol))

#read pkt206
pm1<-read.table(file="PKT206LSSAstate.txt",header=TRUE, sep="\t")
head(pm1)
pm2<-unique(as.character(pm1$node))

common206<-intersect(rd,pm2)
colnames(pm1)=c("V1","V2","V3","V4","V5")
# the index of mp column is symbol,p53wt DD ON, p53wt DD OFF, p53 mu DD ON, p53 mu DD OFF
mp<-pm1[,c(1,2,3)]

mp1<-mp[mp$V1%in%common206,]
write.table(common206,"PKT206wt_mapped_gene.txt",sep="\t",quote=FALSE)
mp1[which((mp1$V2==1)&(mp1$V3==1)),4]=0
mp1[which((mp1$V2==0)&(mp1$V3==0)),4]=0
mp1[which((mp1$V2=="NaN")&(mp1$V3=="NaN")),4]=0
mp1[which((mp1$V2==0)&(mp1$V3==1)),4]=-1
mp1[which((mp1$V2==0)&(mp1$V3=="NaN")),4]=-1
mp1[which((mp1$V2=="NaN")&(mp1$V3==1)),4]=-1
mp1[which((mp1$V2==1)&(mp1$V3==0)),4]=1
mp1[which((mp1$V2==1)&(mp1$V3=="NaN")),4]=1
mp1[which((mp1$V2=="NaN")&(mp1$V3==0)),4]=1
colnames(mp1)<-c("symbol","p53wtDDON","p53wtDDOFF","wt_tr_vs_wt_untr")
x=mp1$wt_tr_vs_wt_untr
table(x)

#identify DE genes from RNA-seq results
d5<-d2[d2$symbol%in%common206,]
d6<-d5[,c(2,5,7)]
d6[,4]=0
d6[which((as.numeric(d6$logFC)>log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=1
d6[which((as.numeric(d6$logFC)< -log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=-1
#d6[which((as.numeric(d6$logFC)>1.5)&(as.numeric(d6$PValue)<0.05)),4]=1
#d6[which((as.numeric(d6$logFC)< -1.5)&(as.numeric(d6$PValue)<0.05)),4]=-1
d6<-d6[!duplicated(d6$symbol),]


d7<-d6[,c(3,4)]
head(mp1)
mp2<-mp1[,c(1,4)]

table(d7$V4)
d8<-merge(d7,mp2, by.x="symbol",by.y="symbol")

colnames(d8)<-c("symbol","Eexp","Emod")

d8$P<-abs(as.numeric(d8$Emod)-as.numeric(d8$Eexp))
x<-d8$P
table(x)


write.table(d8,file="p53wt_tr_vs_wt_untr_206_results",sep="\t",quote=FALSE,row.names = FALSE)

write.table(d7,file="p53wt_tr_vs_wt_untr_Exp_206_results.txt",sep="\t",quote=FALSE,row.names = TRUE)
write.table(mp2,file="p53wt_tr_vs_wt_untr_Mod_206_results.txt",sep="\t",quote=FALSE,row.names = FALSE)

#table(d7$V4)

#p53 mu tr vs p53 wt tr
setwd("D:/work/pkt206")
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)
library("BiocParallel")
library(edgeR)

wt1<-read.table(file="countTableOfp53_wt_tr.txt",header=TRUE, sep="\t")

m2<-read.table(file="countTable_p53muTr.txt",header=TRUE, sep="\t")
#cbind(t, z[, "symbol"][match(rownames(t), rownames(z))])
#countdata1<-cbind(wt1,m2)

countdata1<-cbind(wt1,m2[,"X649PT.bam"][match(rownames(wt1),rownames(m2))])

z<-countdata1[is.na(countdata1),]

group1<-factor(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2))
library(edgeR)
y1<-DGEList(counts=countdata1,group=group1)
#perform the TMM normalization
y1<-calcNormFactors(y1)
design1<-model.matrix(~group1)

y1<-estimateDisp(y1,design1)
bcv<-0.4
#bcv<-0.2
et1<-exactTest(y1,dispersion=bcv^2)
result_step_1_2<-topTags(et1,n=nrow(et1),adjust.method="BH", sort.by="logFC")
#result_step_1_3<-topTags(et1,n=57906,adjust.method="BH", sort.by="logFC")

#annotation
library("org.Hs.eg.db")
d0<-result_step_1_2$table
d1<-cbind(rownames(d0),d0)
colnames(d1)<-c("Id","logFC","logCPM","PValue","FDR")

d1$symbol <- mapIds(org.Hs.eg.db,
                    keys=as.character(d1$Id),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
d2<-na.omit(d1)

write.csv(d2,file="p53_mu_tr_vs_wt_tr.csv",sep="\t",quote=FALSE)
rd<-unique(as.character(d2$symbol))
length(rd)
#read PKT206
pm1<-read.table(file="PKT206LSSAstate.txt",header=TRUE, sep="\t")
head(pm1)
colnames(pm1)=c("V1","V2","V3","V4","V5")
pm2<-unique(as.character(pm1$V1))

common206<-intersect(rd,pm2)
# the index of mp column is symbol, p53 wt DD ON, p53 wt DD OFF, p53 mu DD ON, p53 mu DD OFF
mp<-pm1[,c(1,2,4)]

mp1<-mp[mp$V1%in%common206,]
mp1[which((mp1$V2==1)&(mp1$V4==1)),4]=0
mp1[which((mp1$V2==0)&(mp1$V4==0)),4]=0
mp1[which((mp1$V2=="NaN")&(mp1$V4=="NaN")),4]=0
mp1[which((mp1$V2==0)&(mp1$V4==1)),4]=1
mp1[which((mp1$V2==0)&(mp1$V4=="NaN")),4]=1
mp1[which((mp1$V2=="NaN")&(mp1$V4==1)),4]=1
mp1[which((mp1$V2==1)&(mp1$V4==0)),4]=-1
mp1[which((mp1$V2==1)&(mp1$V4=="NaN")),4]=-1
mp1[which((mp1$V2=="NaN")&(mp1$V4==0)),4]=-1

colnames(mp1)<-c("symbol","p53wtDDON","p53muDDON","mu_tr_vs_wt_tr")
x=mp1$mu_tr_vs_wt_tr
table(x)
#length(common346)
#identify DE genes from RNA-seq results
d5<-d2[d2$symbol%in%common206,]
d6<-d5[,c(2,4,6)]

d6[,4]=0
#d6[which(as.numeric(d6$logFC)>1.5),4]=1
#d6[which(as.numeric(d6$logFC)< -1.5),4]=-1
d6[which((as.numeric(d6$logFC)>log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=1
d6[which((as.numeric(d6$logFC)< -log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=-1
d6=d6[!duplicated(d6$symbol),]

d7<-d6[,c(3,4)]
head(mp1)
mp2<-mp1[,c(1,4)]
table(d7$V4)

d8<-merge(d7,mp2, by.x="symbol",by.y="symbol")

colnames(d8)<-c("symbol","Eexp","Emod")

d8$P<-abs(as.numeric(d8$Emod)-as.numeric(d8$Eexp))
x<-d8$P
table(x)

#c<-d8[duplicated(d8$symbol),]

#write.table(d8,file="p53mu_tr_vs_mu_untr_346_test",sep="\t",quote=FALSE)
write.table(d8,file="p53mu_tr_vs_wt_tr_206_results.txt",sep="\t",quote=FALSE,row.names=FALSE)

write.table(d7,file="p53mu_tr_vs_wt_tr_Exp_206_results.txt",sep="\t",quote=FALSE,row.names = TRUE)
write.table(mp2,file="p53mu_tr_vs_wt_tr_Mod_206_results.txt",sep="\t",quote=FALSE,row.names = FALSE)

##p53 mu tr vs p53 mu untr
setwd("D:/work/pkt206")
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)
library("BiocParallel")
library(edgeR)

m1<-read.table(file="countTableOfp53_mu_untr.txt",header=TRUE, sep="\t")
m2<-read.table(file="countTable_p53muTr.txt",header=TRUE, sep="\t")
m3<-cbind(m1,m2)
colnames(m3)<-c("602PT","667PT","M17PT","M20PT","M37PT","M43PT","M49PT","M50PT","M57PT","M60PT","M61PT","M62PT","M69PT","M77PT","M8PT","M97PT","M99PT","649PT")

countdata1<-m3
group1<-factor(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2))


library(edgeR)
y1<-DGEList(counts=countdata1,group=group1)
#perform the TMM normalization
y1<-calcNormFactors(y1)
design1<-model.matrix(~group1)

y1<-estimateDisp(y1,design1)
bcv<-0.4
#bcv<-0.2
et1<-exactTest(y1,dispersion=bcv^2)
result_step_1_2<-topTags(et1,n=nrow(et1),adjust.method="BH", sort.by="logFC")
#t<-result_step_1_2$table

library("org.Hs.eg.db")
mt0<-result_step_1_2$table
mt1<-cbind(rownames(mt0),mt0)
colnames(mt1)<-c("Id","logFC","logCPM","PValue","FDR")

mt1$symbol <- mapIds(org.Hs.eg.db,
                     keys=as.character(mt1$Id),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
mt2<-na.omit(mt1)
write.csv(mt2,file="p53_mu_tr_vs_mu_untr.csv",sep="\t",quote=FALSE)
rd<-unique(as.character(mt2$symbol))
length(rd)
#read MPH346
pm1<-read.table(file="PKT206LSSAstate.txt",header=TRUE, sep="\t")
head(pm1)
colnames(pm1)=c("V1","V2","V3","V4","V5")
pm2<-unique(as.character(pm1$V1))

common206<-intersect(rd,pm2)
# the index of mp column is symbol, p53 mu DD ON, p53 mu DD OFF
mp<-pm1[,c(1,4,5)]

mp1<-mp[mp$V1%in%common206,]
mp1[which((mp1$V4==1)&(mp1$V5==1)),4]=0
mp1[which((mp1$V4==0)&(mp1$V5==0)),4]=0
mp1[which((mp1$V4=="NaN")&(mp1$V5=="NaN")),4]=0
mp1[which((mp1$V4==0)&(mp1$V5==1)),4]=-1
mp1[which((mp1$V4==0)&(mp1$V5=="NaN")),4]=-1
mp1[which((mp1$V4=="NaN")&(mp1$V5==1)),4]=-1
mp1[which((mp1$V4==1)&(mp1$V5==0)),4]=1
mp1[which((mp1$V4==1)&(mp1$V5=="NaN")),4]=1
mp1[which((mp1$V4=="NaN")&(mp1$V5==0)),4]=1
colnames(mp1)<-c("symbol","p53muDDON","p53muDDOFF","mu_tr_vs_mu_untr")
x=mp1$mu_tr_vs_mu_untr
table(x)
#length(common346)
#identify DE genes from RNA-seq results
d5<-mt2[mt2$symbol%in%common206,]
d6<-d5[,c(2,4,6)]
d6[,4]=0
#d6[which(as.numeric(d6$logFC)>log2(1.5)),4]=1
#d6[which(as.numeric(d6$logFC)< -log2(1.5)),4]=-1

d6[which((as.numeric(d6$logFC)>log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=1
d6[which((as.numeric(d6$logFC)< -log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=-1

d6=d6[!duplicated(d6$symbol),]

d7<-d6[,c(3,4)]
head(mp1)
mp2<-mp1[,c(1,4)]
table(d7$V4)

d8<-merge(d7,mp2, by.x="symbol",by.y="symbol")

colnames(d8)<-c("symbol","Eexp","Emod")

d8$P<-abs(as.numeric(d8$Emod)-as.numeric(d8$Eexp))
x<-d8$P
table(x)

#c<-d8[duplicated(d8$symbol),]

#write.table(d8,file="p53mu_tr_vs_mu_untr_346_test",sep="\t",quote=FALSE)
write.table(d8,file="p53mu_tr_vs_mu_untr_206_results.txt",sep="\t",quote=FALSE,row.names=FALSE)

write.table(d7,file="p53mu_tr_vs_mu_untr_Exp_206_results.txt",sep="\t",quote=FALSE,row.names = TRUE)
write.table(mp2,file="p53mu_tr_vs_mu_untr_Mod_206_results.txt",sep="\t",quote=FALSE,row.names = FALSE)

##p53 wt untr vs p53 mu untr

setwd("D:/work/pkt206")
library(GenomicFeatures)
library(GenomicRanges)
library(Rsamtools)
library(DESeq2)
library(GenomicAlignments)
library("BiocParallel")
library(edgeR)
wunt1<-read.table(file="countTableOfp53_wt_untr.txt",header=TRUE, sep="\t")
m1<-read.table(file="countTableOfp53_mu_untr.txt",header=TRUE, sep="\t")
countdata<-cbind(wunt1,m1)
group1<-factor(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2))

library(edgeR)
y<-DGEList(counts=countdata,group=group1)
#perform the TMM normalization
y<-calcNormFactors(y)
design<-model.matrix(~group1)
#To estimate common dispersion and tagwise dispersions in one run
y<-estimateDisp(y,design)
et<-exactTest(y)
topTags(et)
#use GLM  likelihood ratio test to determing Differential expressed genes
fit<-glmFit(y,design)
# compare mu vs wt
#lrt.2vs1 <- glmLRT(fit, coef=2)
tr <- glmTreat(fit, coef=2)

result_step_1<-topTags(tr,n=nrow(tr$table),adjust.method="BH", sort.by="logFC")
#result_step_2<-topTags(tr,n=57906,adjust.method="BH", sort.by="logFC")


#annotation
library("org.Hs.eg.db")
d0<-result_step_1$table
d1<-cbind(rownames(d0),d0)
colnames(d1)<-c("Id","logFC","unshrunk.logFC","logCPM","PValue","FDR")

d1$symbol <- mapIds(org.Hs.eg.db,
                    keys=as.character(d1$Id),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")


d2<-na.omit(d1)

write.csv(d2,file="p53_wt_untr_vs_p53_mu_untr.csv",sep="\t",quote=FALSE)

rd<-unique(as.character(d2$symbol))

#read MPH346
pm1<-read.table(file="PKT206LSSAstate.txt",header=TRUE, sep="\t")
head(pm1)
colnames(pm1)=c("V1","V2","V3","V4","V5")
pm2<-unique(as.character(pm1$V1))

common206<-intersect(rd,pm2)
# the index of mp column is symbol,p53wt DD ON, p53wt DD OFF, p53 mu DD ON, p53 mu DD OFF
mp<-pm1[,c(1,3,5)]

mp1<-mp[mp$V1%in%common206,]
#write.table(common346,"test_wt_mapped_gene.txt",sep="\t",quote=FALSE)
mp1[which((mp1$V3==1)&(mp1$V5==1)),4]=0
mp1[which((mp1$V3==0)&(mp1$V5==0)),4]=0
mp1[which((mp1$V3=="NaN")&(mp1$V5=="NaN")),4]=0
mp1[which((mp1$V3==0)&(mp1$V5==1)),4]=1
mp1[which((mp1$V3==0)&(mp1$V5=="NaN")),4]=1
mp1[which((mp1$V3=="NaN")&(mp1$V5==1)),4]=1
mp1[which((mp1$V3==1)&(mp1$V5==0)),4]=-1
mp1[which((mp1$V3==1)&(mp1$V5=="NaN")),4]=-1
mp1[which((mp1$V3=="NaN")&(mp1$V5==0)),4]=-1


colnames(mp1)<-c("symbol","p53wtDDOFF","p53muDDOFF","mu_untr_vs_wt_untr")
x=mp1$mu_untr_vs_wt_untr
table(x)

#identify DE genes from RNA-seq results
d5<-d2[d2$symbol%in%common206,]
d6<-d5[,c(2,5,7)]
d6[,4]=0
#d6[which((as.numeric(d6$logFC)>1.5)&(as.numeric(d6$PValue)<0.05)),4]=1
#d6[which((as.numeric(d6$logFC)< -1.5)&(as.numeric(d6$PValue)<0.05)),4]=-1
d6[which((as.numeric(d6$logFC)>log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=1
d6[which((as.numeric(d6$logFC)< -log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=-1
d6<-d6[!duplicated(d6$symbol),]


d7<-d6[,c(3,4)]
head(mp1)
mp2<-mp1[,c(1,4)]

table(d7$V4)
d8<-merge(d7,mp2, by.x="symbol",by.y="symbol")

colnames(d8)<-c("symbol","Eexp","Emod")

d8$P<-abs(as.numeric(d8$Emod)-as.numeric(d8$Eexp))
x<-d8$P
table(x)


write.table(d8,file="p53mu_untr_vs_wt_untr_206_results.txt",sep="\t",quote=FALSE,row.names = FALSE)

write.table(d7,file="p53mu_untr_vs_wt_untr_Exp_206_results.txt",sep="\t",quote=FALSE,row.names = TRUE)
write.table(mp2,file="p53mu_untr_vs_wt_untr_Mod_206_results.txt",sep="\t",quote=FALSE,row.names = FALSE)


#mero14
setwd("D:/work/pkt206/mero14")
library(stringr)
data0=read.csv(file="EtopsideMero14.csv",header=TRUE,sep=",",quote="\"")
d1=data0
rd=unique(as.character(d1$Gene.Symbol))

pm1<-read.table(file="PKT206LSSAstate.txt",header=TRUE, sep="\t")
head(pm1)
colnames(pm1)=c("V1","V2","V3","V4","V5")
pm2<-unique(as.character(pm1$V1))

t<-intersect(rd,pm2)

data1=data0=read.csv(file="GemcitabineMero14.csv",header=TRUE,sep=",",quote="\"")
d2=data1
rd1=unique(as.character(d2$Gene.Symbol))
t2<-intersect(rd1,pm2)
#common206<-intersect(rd,pm2)



mp<-pm1[,c(1,4,5)]

mp1<-mp[mp$V1%in%common206,]
mp1[which((mp1$V4==1)&(mp1$V5==1)),4]=0
mp1[which((mp1$V4==0)&(mp1$V5==0)),4]=0
mp1[which((mp1$V4=="NaN")&(mp1$V5=="NaN")),4]=0
mp1[which((mp1$V4==0)&(mp1$V5==1)),4]=-1
mp1[which((mp1$V4==0)&(mp1$V5=="NaN")),4]=-1
mp1[which((mp1$V4=="NaN")&(mp1$V5==1)),4]=-1
mp1[which((mp1$V4==1)&(mp1$V5==0)),4]=1
mp1[which((mp1$V4==1)&(mp1$V5=="NaN")),4]=1
mp1[which((mp1$V4=="NaN")&(mp1$V5==0)),4]=1
colnames(mp1)<-c("symbol","p53muDDON","p53muDDOFF","mu_tr_vs_mu_untr")
x=mp1$mu_tr_vs_mu_untr
table(x)
#length(common346)
#identify DE genes from RNA-seq results
d5<-mt2[mt2$symbol%in%common206,]
d6<-d5[,c(2,4,6)]
d6[,4]=0
#d6[which(as.numeric(d6$logFC)>log2(1.5)),4]=1
#d6[which(as.numeric(d6$logFC)< -log2(1.5)),4]=-1

d6[which((as.numeric(d6$logFC)>log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=1
d6[which((as.numeric(d6$logFC)< -log2(1.5))&(as.numeric(d6$PValue)<0.05)),4]=-1

d6=d6[!duplicated(d6$symbol),]

d7<-d6[,c(3,4)]
head(mp1)
mp2<-mp1[,c(1,4)]
table(d7$V4)

d8<-merge(d7,mp2, by.x="symbol",by.y="symbol")

colnames(d8)<-c("symbol","Eexp","Emod")

d8$P<-abs(as.numeric(d8$Emod)-as.numeric(d8$Eexp))
x<-d8$P
table(x)



