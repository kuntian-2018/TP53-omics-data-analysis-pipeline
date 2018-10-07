setwd("D:/work/pkt206/mero14_2")

mero14_normalized_input <- read.csv("D:/work/pkt206/mero14_2/mero14_normalized_input.csv")
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
d6$gene_symbol=rownames(d6)
d6_DEgene=d6[which(((as.numeric(d6$logFC)>log2(1.5))&(as.numeric(d6$P.Value)<0.05))|((as.numeric(d6$logFC)< -log2(1.5))&(as.numeric(d6$P.Value)<0.05))),]
#write.csv(d6_DEgene,file="Gem_vs_control_DE_genes.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)
# get up-regulated genes
d6_up=d6[which((as.numeric(d6$logFC)>log2(1.5))&(as.numeric(d6$P.Value)<0.05)),]
#write.csv(d6_up,file="Gem_vs_control_up_regulated_genes.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)
# get down-regulated genes
d6_down=d6[which((as.numeric(d6$logFC)< -log2(1.5))&(as.numeric(d6$P.Value)<0.05)),]
#write.csv(d6_down,file="Gem_vs_control_down_regulated_genes.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)

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
de6$gene_symbol=rownames(de6)
de6_DEgene=de6[which(((as.numeric(de6$logFC)>log2(1.5))&(as.numeric(de6$P.Value)<0.05))|((as.numeric(de6$logFC)< -log2(1.5))&(as.numeric(de6$P.Value)<0.05))),]
#write.csv(de6_DEgene,file="Etop_vs_control_DE_genes.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)
# get up-regulated genes
de6_up=de6[which((as.numeric(de6$logFC)>log2(1.5))&(as.numeric(de6$P.Value)<0.05)),]
#write.csv(de6_up,file="Etop_vs_control_up_regulated_genes.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)
# get down-regulated genes
de6_down=de6[which((as.numeric(de6$logFC)< -log2(1.5))&(as.numeric(de6$P.Value)<0.05)),]
#write.csv(de6_down,file="Etop_vs_control_down_regulated_genes.csv",quote=TRUE,sep="",row.names = FALSE,col.names=TRUE)
#plot Venn diagram
library(VennDiagram)
data1=as.character(d6_DEgene$gene_symbol)
data2=as.character(de6_DEgene$gene_symbol)



data3=as.character(d6_up$gene_symbol)
data4=as.character(de6_up$gene_symbol)
data5=as.character(d6_down$gene_symbol)
data6=as.character(de6_down$gene_symbol)

venn.plot<-draw.quad.venn(
  area1 = length(data3),
  area2 = length(data4), 
  area3 = length(data5),
  area4 = length(data6),
  n12=length(intersect(data3,data4)), 
  n13=length(intersect(data3,data5)),
  n14=length(intersect(data3,data6)),
  n23=length(intersect(data4,data5)),
  n24=length(intersect(data4,data6)),
  n34=length(intersect(data5,data6)),
  n123=length(Reduce(intersect, list(data3,data4,data5))),
  n124=length(Reduce(intersect, list(data3,data4,data6))),
  n134=length(Reduce(intersect, list(data3,data5,data6))),
  n234=length(Reduce(intersect, list(data4,data5,data6))),
  n1234=length(Reduce(intersect, list(data3,data4,data5,data6))),
  category = c("Gem_up", "Etop_up","Gem_down","Etop_down"),
  fill = c("orange","red", "green","blue"),
  lty="solid",
  cex =2, 
  cat.cex=2,
  cat.col =c("orange","red", "green","blue")
)


common_gem_down_and_etop_down=intersect(data5,data6)
write.table(common_gem_down_and_etop_down,file="Gem_down_vs_Etop_down_regulated_overlap.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

common_gem_up_and_etop_up=intersect(data3,data4)
write.table(common_gem_up_and_etop_up,file="Gem_up_vs_Etop_up_regulated_overlap.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

common_gem_down_and_etop_up=intersect(data5,data4)

uniq_gem_up=setdiff(data3,common_gem_up_and_etop_up)

write.table(uniq_gem_up,file="Only_up_regulated_in_Gem.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

uniq_gem_down=setdiff(data5, common_gem_down_and_etop_down)
write.table(uniq_gem_down,file="Only_down_regulated_in_Gem.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


uniq_etop_down=setdiff(data6,common_gem_down_and_etop_down)

write.table(uniq_etop_down,file="Only_down_regulated_in_Etop.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

common_a1=union(common_gem_up_and_etop_up,common_gem_down_and_etop_up)
uniq_etop_up=setdiff(data4,common_a1)
write.table(uniq_etop_up,file="Only_up_regulated_in_Etop.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



#write.table(common_gem_down_and_etop_up,file="Gem_down_vs_Etop_up_regulated_overlap.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

#df1<-setdiff(data1,data2)
#df2<-setdiff(data2,data1)
#write.table(df1,file="Unique_gene_in_Gem.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
#write.table(df2,file="Unique_gene_in_Etop.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

#out.plot <-  draw.pairwise.venn(
#  area1 = length(data1),
#  area2 = length(data2),
#  cross.area = length(intersect(data1,data2)),
#  category = c("Gem_DE_gene", "Etop_DE_gene"),
#  fill = c("orange", "red"),
#  lty = "blank",
#  cex = 2,
#  cat.cex = 1.3,
#  cat.pos = c(0, 0),
#  cat.dist = 0.03,
#  cat.just = list(c(0, 0), c(0, 0)),
#  ext.pos = 30,
#  ext.dist = -0.05,
#  ext.length = 0.85,
#  ext.line.lwd = 2,
#  ext.line.lty = "dashed"
#)
#common_RMA=intersect(data1,data2)
#write.table(common_RMA,file="Gem_vs_Etop_RMA_overlap.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

#df1<-setdiff(data1,data2)
#df2<-setdiff(data2,data1)
#write.table(df1,file="Unique_gene_in_Gem_vs_control.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
#write.table(df2,file="Unique_gene_in_Etop_vs_control.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)






