setwd("D:/work/pkt206/mero14_3")
mero14_normalized_input <- read.csv("D:/work/pkt206/mero14_3/mero14_normalized_input.csv")
input1=mero14_normalized_input
head(input1)
#retrieve columns of probe id,control signal, gem signal, etopside signal and gene symbol
input2=input1[,c(1,3,5,6,8,10,11,14)]
head(input2)
input3=input2[,c(8,2:7)]
#write.csv(input3,file="retrieved_columns.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)
nrow(input3)
## remove rows with empty gene symbol value,"---"
input4=input3[!(input3$Gene.Symbol=="---"),]
nrow(input4)
#write.csv(input4,file="empty_gene_symbol_removed.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)
##calculate median value of genes for each sample
d1=input4
d2=d1
d3=aggregate(d2[,2:7],list(d2$Gene.Symbol),median)
#length(unique(d1$Gene.Symbol))
#write.csv(d3,file="median_value_for_same_gene.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)
#control group: A2824.05,A2824.10;which are in column 4 and column 7
#gem group: A2824.02, A2824.07; which are in column 2 and column 5
#Etopside group:A2824.04, A2824.09,which are in column 3 and column 6.

#prepare expression matrix of Gem and control
d4=as.matrix(d3[,c(4,7,2,5)])
rownames(d4)=d3$Group.1
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biobase")
library("Biobase")
expr <- ExpressionSet(assayData=d4)
library(limma)
design2 <- model.matrix(~0+factor(c(1,1,2,2)))
colnames(design2) <- c("Control", "Gem")
contrast.matrix2 <- makeContrasts(Gem-Control, levels=design2)
fit <- lmFit(expr,design2)

fit2 <- contrasts.fit(fit, contrast.matrix2)

ebfit2 <- eBayes(fit2)
# log2-fold change threshold of 1.5
result <- topTable(ebfit2, coef=1, number=nrow(expr))
#write.table(result,file="topTable_of_Gem_vs_control_mero14.txt",sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)
# filter DE genes by cut off of Fold change over 1.5 or -1.5 and p value less than 0.05
d6=result
d7=d6[,c(1,4)]
d7$gene_symbol=rownames(d7)

STSFA_input <- read.csv("D:/work/pkt206/mero14_3/STSFA_input.csv")
m1=STSFA_input
m1$gem_FC=(m1$medium_gem/m1$medium_control)
m1$etop_FC=(m1$meadium_etop/m1$medium_control)
m2=m1[,c(1,5,6)]

#write.csv(m2,file="STSFA_mero14_FC_stage_1.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
m3=na.omit(m2)
m3$gem_lfc=log10(m3$gem_FC)
m3$etop_lfc=log10(m3$etop_FC)

#write.csv(m3,file="STSFA_mero14_FC_stage_2.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
gem_mean=mean(m3$gem_lfc)
gem_std=sd(m3$gem_lfc)
etop_mean=mean(m3$etop_lfc)
etop_std=sd(m3$etop_lfc)
gem_up=gem_mean+gem_std
gem_low=gem_mean-gem_std
etop_up=etop_mean+etop_std
etop_low=etop_mean-etop_std
m3$gem_mod=0
m3$etop_mod=0

m3[m3$gem_lfc>gem_up,6]=1

m3[m3$gem_lfc<gem_low,6]=-1
m3[m3$etop_lfc>etop_up ,7]=1

m3[m3$etop_lfc< etop_low,7]=-1
table(m3$gem_mod)
#-1,0,1
#20, 157,19
table(m3$etop_mod)
#-1,0,1
#15,163,18

write.csv(m3,file="STSFA_mero14_FC_stage_3.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
m4=m3[,c(1,6,7)]
#write.csv(m4,file="STSFA_mero14_state_value.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)


#compare STSFA with microarry results 

# Gem_vs_control
common_gene=intersect(m4$Gene,d7$gene_symbol)
d8=d7[which(d7$gene_symbol%in%common_gene),]
d8$exp=0

d8[which((as.numeric(d8$logFC)>log2(1.5))&(as.numeric(d8$P.Value)<0.05)),4]=1
d8[which((as.numeric(d8$logFC)< -log2(1.5))&(as.numeric(d8$P.Value)<0.05)),4]=-1
#write.csv(d8,file="gem_vs_control_mero14_exp.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
d9<-d8[,c(3,4)]

mp2<-m3[,c(1,6)]

table(d9$exp)
# -1,0,1
#12,152,27
d10<-merge(d9,mp2, by.x="gene_symbol",by.y="Gene")

colnames(d10)<-c("symbol","Eexp","Emod")

d10$P<-abs(as.numeric(d10$Emod)-as.numeric(d10$Eexp))
x<-d10$P
table(x)
#0,1
#163,28
#write.csv(d10,file="STSFA_vs_microarray_of_gem_vs_control.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)


#prepare expression matrix of Etop and control
de1=as.matrix(d3[,c(4,7,3,6)])
rownames(de1)=d3$Group.1
expr1 <- ExpressionSet(assayData=de1)
#library(limma)
design3 <- model.matrix(~0+factor(c(1,1,2,2)))
colnames(design3) <- c("Control", "Etop")
contrast.matrix3 <- makeContrasts(Etop-Control, levels=design3)
fit3 <- lmFit(expr1,design3)

fit4 <- contrasts.fit(fit3, contrast.matrix3)

ebfit4 <- eBayes(fit4)
# log2-fold change threshold of 1.5
result_2 <- topTable(ebfit4, coef=1, number=nrow(expr1))
#write.table(result_2,file="topTable_of_Etop_vs_control_mero14.txt",sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)
de6=result_2

de7=de6[,c(1,4)]
de7$gene_symbol=rownames(de7)
common_gene_1=intersect(m4$Gene,de7$gene_symbol)


de8=de7[which(de7$gene_symbol%in%common_gene_1),]
de8$exp=0

de8[which((as.numeric(de8$logFC)>log2(1.5))&(as.numeric(de8$P.Value)<0.05)),4]=1
de8[which((as.numeric(de8$logFC)< -log2(1.5))&(as.numeric(de8$P.Value)<0.05)),4]=-1


write.csv(de8,file="etop_vs_control_mero14_exp.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
de9<-de8[,c(3,4)]

mp3<-m3[,c(1,7)]

table(de9$exp)
# -1,0,1
#14,146,31
de10<-merge(de9,mp3, by.x="gene_symbol",by.y="Gene")

colnames(de10)<-c("symbol","Eexp","Emod")

de10$P<-abs(as.numeric(de10$Emod)-as.numeric(de10$Eexp))
xe<-de10$P
table(xe)
#0,1
#157,34
write.csv(de10,file="STSFA_vs_microarray_of_etop_vs_control.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)


