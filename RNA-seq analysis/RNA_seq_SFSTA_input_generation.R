setwd("D:/work/pkt206/sfsta")
library(stringr)


#m1<-read.table(file="wtPKT206.txt",header=FALSE, sep="")
m1<-read.table(file="muPKT206.txt",header=FALSE, sep="")
c1<-m1
c4<-c1
c5=c4[order(c4[,1]),]
#write.table(c5,file="PKT206_ordered.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
write.table(c5,file="muPKT206_ordered.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
#df[!(duplicated(df[c("c","d")]) | duplicated(df[c("c","d")], fromLast = TRUE)), ]
c6=c5[!(duplicated(c5[c(1,3)]) | duplicated(c5[c(1,3)],fromLast=TRUE)),]

c7=c5[(duplicated(c5[c(1,3)]) | duplicated(c5[c(1,3)],fromLast=TRUE)),]

#write.table(c6,file="PKT206_no_ambivalent_edges.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

#write.table(c7,file="PKT206_ambivalent_edges.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
write.table(c6,file="muPKT206_no_ambivalent_edges.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

write.table(c7,file="muPKT206_ambivalent_edges.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
#c6 includes all no ambivalent edges
c8=c6
c8[c8[,2]=="Inhibits",4]=c("inhibition")
c8[c8[,2]=="Activates",4]=c("activation")
c8[,5]=0
c9=c8[,c(1,4,3,5)]
colnames(c9)=c("Nodes1","Type","Nodes2","Weight")

#write.csv(c9,file="wt_PKT206_network.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
write.csv(c9,file="mu_PKT206_network.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)

node0=unique(c(as.character(c9[,1]),as.character(c9[,3])))
write.table(node0,file="muPKT206_node_left.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

m2<-read.table(file="wtPKT206.txt",header=FALSE, sep="")
cw4<-m2
cw5=cw4[order(cw4[,1]),]
cw6=cw5[!(duplicated(cw5[c(1,3)]) | duplicated(cw5[c(1,3)],fromLast=TRUE)),]
node1=unique(c(as.character(cw6[,1]),as.character(cw6[,3])))
#node1 is 206 nodes
node_mu_isolated=setdiff(node1,node0)
write.table(node1,file="wtPKT206_node.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
write.table(node_mu_isolated,file="muPKT206_node_isolated.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

n1=cbind(node1,node1)
n1=as.data.frame(n1)
n1[,3]=c("gene")
n1[,4]=c("no")
n1[,5]=0
colnames(n1)=c("KEY","ENTREZ_ID","NODE_TYPE","TARGET_PROCESS","SCORE")
write.csv(n1,file="wt_PKT206_node_attributes.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
#prepare RNA-seq input file

#wunt1<-read.table(file="countTableOfp53_wt_untr.txt",header=TRUE, sep="\t")
#wtr1<-read.table(file="countTableOfp53_wt_tr.txt",header=TRUE, sep="\t")
mt1<-read.table(file="countTableOfp53_mu_untr.txt",header=TRUE, sep="\t")
mt2<-read.table(file="countTable_p53muTr.txt",header=TRUE, sep="\t")
#mt=cbind(mt1,mt2)
countdata<-cbind(mt1,mt2[,"X649PT.bam"][match(rownames(mt1),rownames(mt2))])
write.csv(countdata,file="combindMuCountTable.csv",sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)

temp<-countdata+2


# log transform the raw read count 
nor_countdata<-log2(temp)*100
#get the normalized count matrix of p53 mutant untreated


#annotation
library("org.Hs.eg.db")
d0<-nor_countdata
d1<-cbind(rownames(d0),d0)
#colnames(d1)<-c("Id","M47PT","M48PT","M52PT","M53PT","M58PT","M608PT","M63PT","M67PT","M68PT","M691PT","M697PT","M699PT","M6PT","M700PT","M701PT","M704PT","M70PT","M710PT","M71PT","M73PT","M75PT","M76PT","M80PT","M82PT","M94PT","M98PT")
#  colnames(d1)<-c("Id","617PT","618PT","634PT","655PT","M101PT","M614PT","M626PT","M632PT","M637PT","M640PT","M645PT","M663PT","M664PT","M666PT","M668PT","M669PT","M684PT","M685PT","M686PT","M687PT","M689PT","M690PT","M694PT","M703PT","M705PT","M714PT","M719PT")
colnames(d1)<-c("Id","602PT","667PT","M17PT","M20PT","M37PT","M43PT","M49PT","M50PT","M57PT","M60PT","M61PT","M62PT","M69PT","M77PT","M8PT","M97PT","M99PT","649PT")
d1$symbol <- mapIds(org.Hs.eg.db,
                    keys=as.character(d1$Id),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")


d2<-na.omit(d1)
#d3=d2[!duplicated(d2$symbol),]
#d4=d2[duplicated(d2$symbol),]
t3<-as.character(d2$symbol)
write.csv(d2,file="p53_Mu_MappedMuCountTable.csv",sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)

node0=read.table(file="wtPKT206_node.txt",header=FALSE, sep="")
node1=as.character(unique(node0[,1]))
common=intersect(node1,t3)
isolated_node_1=setdiff(node1,t3)
#write.table(isolated_node_1,file="RNAseq_to_muPKT206_unmapped_node.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


d3<-d2[d2$symbol%in%common,]
d4<-d3[duplicated(d3$symbol),]
#write.table(d4,file="RNAseq_duplicated_node_in_PKT206.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
#produce sample input ifle
mu_tr_sample=d3[,c(ncol(d3),ncol(d3),ncol(d3)-1)]

write.csv(mu_tr_sample,file="Mu_tr_sample_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
duplicated_name=mu_tr_sample[duplicated(mu_tr_sample$symbol),]

mu_untr_sample_1=d3[,c(ncol(d3),ncol(d3),2)]

write.csv(mu_untr_sample_1,file="Mu_untr_sample_1_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

mu_untr_sample_2=d3[,c(ncol(d3),ncol(d3),3)]

write.csv(mu_untr_sample_2,file="Mu_untr_sample_2_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


mu_untr_sample_3=d3[,c(ncol(d3),ncol(d3),4)]

write.csv(mu_untr_sample_3,file="Mu_untr_sample_3_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


mu_untr_sample_4=d3[,c(ncol(d3),ncol(d3),5)]

write.csv(mu_untr_sample_4,file="Mu_untr_sample_4_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



mu_untr_sample_5=d3[,c(ncol(d3),ncol(d3),6)]

write.csv(mu_untr_sample_5,file="Mu_untr_sample_5_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


mu_untr_sample_6=d3[,c(ncol(d3),ncol(d3),7)]

write.csv(mu_untr_sample_6,file="Mu_untr_sample_6_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)




mu_untr_sample_7=d3[,c(ncol(d3),ncol(d3),8)]

write.csv(mu_untr_sample_7,file="Mu_untr_sample_7_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



mu_untr_sample_8=d3[,c(ncol(d3),ncol(d3),9)]

write.csv(mu_untr_sample_8,file="Mu_untr_sample_8_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



mu_untr_sample_9=d3[,c(ncol(d3),ncol(d3),10)]

write.csv(mu_untr_sample_9,file="Mu_untr_sample_9_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


mu_untr_sample_10=d3[,c(ncol(d3),ncol(d3),11)]

write.csv(mu_untr_sample_10,file="Mu_untr_sample_10_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

mu_untr_sample_11=d3[,c(ncol(d3),ncol(d3),12)]

write.csv(mu_untr_sample_11,file="Mu_untr_sample_11_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)




mu_untr_sample_12=d3[,c(ncol(d3),ncol(d3),13)]

write.csv(mu_untr_sample_12,file="Mu_untr_sample_12_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


mu_untr_sample_13=d3[,c(ncol(d3),ncol(d3),14)]

write.csv(mu_untr_sample_13,file="Mu_untr_sample_13_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


mu_untr_sample_14=d3[,c(ncol(d3),ncol(d3),15)]

write.csv(mu_untr_sample_14,file="Mu_untr_sample_14_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



mu_untr_sample_15=d3[,c(ncol(d3),ncol(d3),16)]

write.csv(mu_untr_sample_15,file="Mu_untr_sample_15_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


mu_untr_sample_16=d3[,c(ncol(d3),ncol(d3),17)]

write.csv(mu_untr_sample_16,file="Mu_untr_sample_16_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



mu_untr_sample_17=d3[,c(ncol(d3),ncol(d3),18)]

write.csv(mu_untr_sample_17,file="Mu_untr_sample_17_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)




## prepare p53 wt treated sample microarray input

setwd("D:/work/pkt206/sfsta")
library(stringr)

wt1<-read.table(file="countTableOfp53_wt_tr.txt",header=TRUE, sep="\t")
#mt2<-read.table(file="countTable_p53muTr.txt",header=TRUE, sep="\t")
#mt=cbind(mt1,mt2)
countdata<-wt1
#write.csv(countdata,file="P53Wt_combindMuCountTable.csv",sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)

temp<-countdata+2


# log transform the raw read count 
nor_countdata<-log2(temp)*100
#get the normalized count matrix of p53 mutant untreated


#annotation
library("org.Hs.eg.db")
d0<-nor_countdata
d1<-cbind(rownames(d0),d0)
#colnames(d1)<-c("Id","M47PT","M48PT","M52PT","M53PT","M58PT","M608PT","M63PT","M67PT","M68PT","M691PT","M697PT","M699PT","M6PT","M700PT","M701PT","M704PT","M70PT","M710PT","M71PT","M73PT","M75PT","M76PT","M80PT","M82PT","M94PT","M98PT")
#each column is the sample name of p53 wt patient
colnames(d1)<-c("Id","617PT","618PT","634PT","655PT","M101PT","M614PT","M626PT","M632PT","M637PT","M640PT","M645PT","M663PT","M664PT","M666PT","M668PT","M669PT","M684PT","M685PT","M686PT","M687PT","M689PT","M690PT","M694PT","M703PT","M705PT","M714PT","M719PT")
#colnames(d1)<-c("Id","602PT","667PT","M17PT","M20PT","M37PT","M43PT","M49PT","M50PT","M57PT","M60PT","M61PT","M62PT","M69PT","M77PT","M8PT","M97PT","M99PT","649PT")
d1$symbol <- mapIds(org.Hs.eg.db,
                    keys=as.character(d1$Id),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")


d2<-na.omit(d1)
#d3=d2[!duplicated(d2$symbol),]
#d4=d2[duplicated(d2$symbol),]
t3<-as.character(d2$symbol)
head(d2)
write.csv(d2,file="p53_wt_MappedMuCountTable.csv",sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)

node0=read.table(file="wtPKT206_node.txt",header=FALSE, sep="")
node1=as.character(unique(node0[,1]))
common=intersect(node1,t3)
isolated_node_1=setdiff(node1,t3)
write.table(isolated_node_1,file="RNAseq_to_wt_PKT206_unmapped_node.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


d3<-d2[d2$symbol%in%common,]
d4<-d3[duplicated(d3$symbol),]

wt_tr_sample_1=d3[,c(ncol(d3),ncol(d3),2)]

write.csv(wt_tr_sample_1,file="Wt_tr_sample_1_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_2=d3[,c(ncol(d3),ncol(d3),3)]

write.csv(wt_tr_sample_2,file="Wt_tr_sample_2_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



wt_tr_sample_3=d3[,c(ncol(d3),ncol(d3),4)]

write.csv(wt_tr_sample_3,file="Wt_tr_sample_3_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



wt_tr_sample_4=d3[,c(ncol(d3),ncol(d3),5)]

write.csv(wt_tr_sample_4,file="Wt_tr_sample_4_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



wt_tr_sample_5=d3[,c(ncol(d3),ncol(d3),6)]

write.csv(wt_tr_sample_5,file="Wt_tr_sample_5_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



wt_tr_sample_6=d3[,c(ncol(d3),ncol(d3),7)]

write.csv(wt_tr_sample_6,file="Wt_tr_sample_6_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_7=d3[,c(ncol(d3),ncol(d3),8)]

write.csv(wt_tr_sample_7,file="Wt_tr_sample_7_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



wt_tr_sample_8=d3[,c(ncol(d3),ncol(d3),9)]

write.csv(wt_tr_sample_8,file="Wt_tr_sample_8_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



wt_tr_sample_9=d3[,c(ncol(d3),ncol(d3),10)]

write.csv(wt_tr_sample_9,file="Wt_tr_sample_9_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_10=d3[,c(ncol(d3),ncol(d3),11)]

write.csv(wt_tr_sample_10,file="Wt_tr_sample_10_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



wt_tr_sample_11=d3[,c(ncol(d3),ncol(d3),12)]

write.csv(wt_tr_sample_11,file="Wt_tr_sample_11_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_12=d3[,c(ncol(d3),ncol(d3),13)]

write.csv(wt_tr_sample_12,file="Wt_tr_sample_12_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)



wt_tr_sample_13=d3[,c(ncol(d3),ncol(d3),14)]

write.csv(wt_tr_sample_13,file="Wt_tr_sample_13_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_14=d3[,c(ncol(d3),ncol(d3),15)]

write.csv(wt_tr_sample_14,file="Wt_tr_sample_14_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_15=d3[,c(ncol(d3),ncol(d3),16)]

write.csv(wt_tr_sample_15,file="Wt_tr_sample_15_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_16=d3[,c(ncol(d3),ncol(d3),17)]

write.csv(wt_tr_sample_16,file="Wt_tr_sample_16_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_17=d3[,c(ncol(d3),ncol(d3),18)]

write.csv(wt_tr_sample_17,file="Wt_tr_sample_17_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_18=d3[,c(ncol(d3),ncol(d3),19)]

write.csv(wt_tr_sample_18,file="Wt_tr_sample_18_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_19=d3[,c(ncol(d3),ncol(d3),20)]

write.csv(wt_tr_sample_19,file="Wt_tr_sample_19_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_20=d3[,c(ncol(d3),ncol(d3),21)]

write.csv(wt_tr_sample_20,file="Wt_tr_sample_20_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_21=d3[,c(ncol(d3),ncol(d3),22)]

write.csv(wt_tr_sample_21,file="Wt_tr_sample_21_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_22=d3[,c(ncol(d3),ncol(d3),23)]

write.csv(wt_tr_sample_22,file="Wt_tr_sample_22_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_23=d3[,c(ncol(d3),ncol(d3),24)]

write.csv(wt_tr_sample_23,file="Wt_tr_sample_23_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_24=d3[,c(ncol(d3),ncol(d3),25)]

write.csv(wt_tr_sample_24,file="Wt_tr_sample_24_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_25=d3[,c(ncol(d3),ncol(d3),26)]

write.csv(wt_tr_sample_25,file="Wt_tr_sample_25_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_26=d3[,c(ncol(d3),ncol(d3),27)]

write.csv(wt_tr_sample_26,file="Wt_tr_sample_26_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_27=d3[,c(ncol(d3),ncol(d3),28)]

write.csv(wt_tr_sample_27,file="Wt_tr_sample_27_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)
## prepare input of wt untreated samples



setwd("D:/work/pkt206/sfsta")
library(stringr)

wt1<-read.table(file="countTableOfp53_wt_untr.txt",header=TRUE, sep="\t")
#mt2<-read.table(file="countTable_p53muTr.txt",header=TRUE, sep="\t")
#mt=cbind(mt1,mt2)
countdata<-wt1
#write.csv(countdata,file="P53Wt_combindMuCountTable.csv",sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)

temp<-countdata+2


# log transform the raw read count 
nor_countdata<-log2(temp)*100
#get the normalized count matrix of p53 mutant untreated


#annotation
library("org.Hs.eg.db")
d0<-nor_countdata
d1<-cbind(rownames(d0),d0)
colnames(d1)<-c("Id","M47PT","M48PT","M52PT","M53PT","M58PT","M608PT","M63PT","M67PT","M68PT","M691PT","M697PT","M699PT","M6PT","M700PT","M701PT","M704PT","M70PT","M710PT","M71PT","M73PT","M75PT","M76PT","M80PT","M82PT","M94PT","M98PT")
#each column is the sample name of p53 wt patient
#colnames(d1)<-c("Id","617PT","618PT","634PT","655PT","M101PT","M614PT","M626PT","M632PT","M637PT","M640PT","M645PT","M663PT","M664PT","M666PT","M668PT","M669PT","M684PT","M685PT","M686PT","M687PT","M689PT","M690PT","M694PT","M703PT","M705PT","M714PT","M719PT")
#colnames(d1)<-c("Id","602PT","667PT","M17PT","M20PT","M37PT","M43PT","M49PT","M50PT","M57PT","M60PT","M61PT","M62PT","M69PT","M77PT","M8PT","M97PT","M99PT","649PT")
d1$symbol <- mapIds(org.Hs.eg.db,
                    keys=as.character(d1$Id),
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")


d2<-na.omit(d1)
#d3=d2[!duplicated(d2$symbol),]
#d4=d2[duplicated(d2$symbol),]
t3<-as.character(d2$symbol)
head(d2)
write.csv(d2,file="p53_wt_untr_MappedCountTable.csv",sep="\t",quote=FALSE,row.names = TRUE,col.names=TRUE)

node0=read.table(file="wtPKT206_node.txt",header=FALSE, sep="")
node1=as.character(unique(node0[,1]))
common=intersect(node1,t3)
isolated_node_1=setdiff(node1,t3)
#write.table(isolated_node_1,file="RNAseq_to_wt_PKT206_unmapped_node.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


d3<-d2[d2$symbol%in%common,]
d4<-d3[duplicated(d3$symbol),]

wt_tr_sample_1=d3[,c(ncol(d3),ncol(d3),2)]

write.csv(wt_tr_sample_1,file="Wt_untr_sample_1_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_2=d3[,c(ncol(d3),ncol(d3),3)]

write.csv(wt_tr_sample_2,file="Wt_untr_sample_2_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_3=d3[,c(ncol(d3),ncol(d3),4)]

write.csv(wt_tr_sample_3,file="Wt_untr_sample_3_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_4=d3[,c(ncol(d3),ncol(d3),5)]

write.csv(wt_tr_sample_4,file="Wt_untr_sample_4_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_5=d3[,c(ncol(d3),ncol(d3),6)]

write.csv(wt_tr_sample_5,file="Wt_untr_sample_5_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_6=d3[,c(ncol(d3),ncol(d3),7)]

write.csv(wt_tr_sample_6,file="Wt_untr_sample_6_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_7=d3[,c(ncol(d3),ncol(d3),8)]

write.csv(wt_tr_sample_7,file="Wt_untr_sample_7_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_8=d3[,c(ncol(d3),ncol(d3),9)]

write.csv(wt_tr_sample_8,file="Wt_untr_sample_8_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_9=d3[,c(ncol(d3),ncol(d3),10)]

write.csv(wt_tr_sample_9,file="Wt_untr_sample_9_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_10=d3[,c(ncol(d3),ncol(d3),11)]

write.csv(wt_tr_sample_10,file="Wt_untr_sample_10_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_11=d3[,c(ncol(d3),ncol(d3),12)]

write.csv(wt_tr_sample_11,file="Wt_untr_sample_11_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_12=d3[,c(ncol(d3),ncol(d3),13)]

write.csv(wt_tr_sample_12,file="Wt_untr_sample_12_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_13=d3[,c(ncol(d3),ncol(d3),14)]

write.csv(wt_tr_sample_13,file="Wt_untr_sample_13_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_14=d3[,c(ncol(d3),ncol(d3),15)]

write.csv(wt_tr_sample_14,file="Wt_untr_sample_14_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_15=d3[,c(ncol(d3),ncol(d3),16)]

write.csv(wt_tr_sample_15,file="Wt_untr_sample_15_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_16=d3[,c(ncol(d3),ncol(d3),17)]

write.csv(wt_tr_sample_16,file="Wt_untr_sample_16_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)


wt_tr_sample_17=d3[,c(ncol(d3),ncol(d3),18)]

write.csv(wt_tr_sample_17,file="Wt_untr_sample_17_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_18=d3[,c(ncol(d3),ncol(d3),19)]

write.csv(wt_tr_sample_18,file="Wt_untr_sample_18_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_19=d3[,c(ncol(d3),ncol(d3),20)]

write.csv(wt_tr_sample_19,file="Wt_untr_sample_19_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_20=d3[,c(ncol(d3),ncol(d3),21)]

write.csv(wt_tr_sample_20,file="Wt_untr_sample_20_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_21=d3[,c(ncol(d3),ncol(d3),22)]

write.csv(wt_tr_sample_21,file="Wt_untr_sample_21_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_22=d3[,c(ncol(d3),ncol(d3),23)]

write.csv(wt_tr_sample_22,file="Wt_untr_sample_22_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_23=d3[,c(ncol(d3),ncol(d3),24)]

write.csv(wt_tr_sample_23,file="Wt_untr_sample_23_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_24=d3[,c(ncol(d3),ncol(d3),25)]

write.csv(wt_tr_sample_24,file="Wt_untr_sample_24_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_25=d3[,c(ncol(d3),ncol(d3),26)]

write.csv(wt_tr_sample_25,file="Wt_untr_sample_25_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)

wt_tr_sample_26=d3[,c(ncol(d3),ncol(d3),27)]

write.csv(wt_tr_sample_26,file="Wt_untr_sample_26_input.csv",sep="\t",quote=FALSE,row.names = FALSE,col.names=FALSE)






