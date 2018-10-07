

#source("https://bioconductor.org/biocLite.R")
#biocLite("RTNsurvival")
#library(RTNsurvival)
install.packages(c("survival", "survminer"))
setwd("D:/work/pkt206/Cox")
library("survival")
library("survminer")
#data=read.csv(file="patient_data_1.csv",header=TRUE)
#data=read.csv(file="patient_data_65_samples.csv",header=TRUE)
data=read.csv(file="patient_data_71_samples.csv",header=TRUE)
#data=read.csv(file="patient_data_71_samples_test_stage.csv",header=TRUE)
data1=data[,-1]
#test=data[duplicated(data$Node_Name),]


rownames(data1)=data$Node_Name
head(data1)
data1$survival=round(data1$survival*365)
data2=data1
#data1[,4:ncol(data2)]=log10(abs(data1[,4:ncol(data1)]+10))
#data1[,4:ncol(data2)]=log2(abs(data1[,4:ncol(data1)]+1))
data1[,4:ncol(data2)]=data1[,4:ncol(data1)]/100
#write.csv(data1,file="testlog10.csv")


#multivariate Cox
#test all 203 genes
#covariates <- c("AATF","ABCB1","ABCC1","AIFM2","APAF1","AR","ARID3A","ATF3","ATM","ATR","AURKA","AXIN1","Apoptosis","BAK1","BAX","BBC3","BCL2","BCL3","BCL6","BDKRB1","BNIP3L",	"BRCA1","BTG2","CALD1","CASP8","CCNA2","CCNB1","CCND1","CCNG1","CD44","CD58","CD59","CD82","CDC20","CDC25A","CDK2","CDK4","CDK5","CDK9","CDKN1A","CDKN1B","CDKN2A","CHEK1","CHEK2","CIAPIN1","CKB","CKM","CKS2","COL18A1","CSNK2A2","CXCR4","Cellularsenescence","DDB2","DDIT4","DDX5","DFNA5","DKK1","DNAdamage","DUSP2","DUSP4","DUSP5","DYRK2","E2F1","ECT2","EDA2R","EGFR","EIF2AK2","ELAVL1","EPHB4","ERBB2","ESR1","EZH2","FAS","FDXR","FEN1","FGF2","FHL2","FOS","FOXM1","GADD45A","GAPDH","GSTP1","GTSE1","H2AFZ","HDAC1","HIC1","HIF1A","HIPK2",	"HIPK4","HMMR","HNF4A","HOXA11","HSP90AB1","HSPA4","HTATIP2","ICAM1","ID3","IER3","IFI16","IFITM2","IFNA1","IGF1R",	"IGFBP1","IGFBP7","IL6","IQCB1","ISG15","KAT2B","KLF4","KRT19","KRT8","LATS2","LTF","MAP4","MAP4K4","MAPK1","MAPK14","MAPK8","MAPK9","MCL1","MCTS1","MDM2",	"MDM4",	"MGMT",	"MMP1",	"MMP13",	"MMP2",	"MSH2",	"MUC1",	"MYC",	"MYCN",	"NCL",	"NLRC4",	"NME1",	"NOTCH1",	"NOV",	"NR2C1",	"NTN1",	"PADI4",	"PARK2",	"PCBP4",	"PCNA",	"PDGFRB",	"PDRG1",	"PEG3",	"PERP",	"PIDD1",	"PLAUR",	"POU4F1",	"PPM1A",	"PPM1D",	"PRC1",	"PRKCA",	"PRKD1",	"PRKDC",	"PRKG1",	"PRSS50",	"PSEN1",	"PSMD10",	"PTEN",	"PTGS2",	"PTTG1",	"RAD51",	"RAF1",	"RAS",	"RECQL4",	"RGCC",	"RGS16",	"RPRM",	"RREB1",	"RRM2B",	"S100A2",	"S100A6",	"SEMA3B",	"SERPINB5",	"SERPINF1",	"SESN2",	"SFN",	"SGK1",	"SIAH1",	"SIVA1",	"SLC2A1",	"SLC2A4",	"SLC6A6",	"SOX4",	"SP7",	"TCF7L2",	"TFDP1",	"TGFA",	"TGFB1",	"THBS1",	"TIAF1",	"TIGAR",	"TLR3",	"TNFRSF10A",	"TNFRSF10B","TP53AIP1",	"TP53I13",	"TP53INP1",	"VEGFA",	"VRK1",	"WWP1",	"XAF1",	"YBX1",	"ZMAT3")
#x=paste(covariates,collapse ='+')


#test multiple variate of all 202 genes
res.cox.test <- coxph(formula=Surv(data1$survival, data1$status) ~ AATF+ABCB1+ABCC1+AIFM2+APAF1+AR+ARID3A+ATF3+ATM+ATR+AURKA+AXIN1+Apoptosis+BAK1+BAX+BBC3+BCL2+BCL3+BCL6+BDKRB1+BNIP3L+BRCA1+BTG2+CALD1+CASP8+CCNA2+CCNB1+CCND1+CCNG1+CD44+CD58+CD59+CD82+CDC20+CDC25A+CDK2+CDK4+CDK5+CDK9+CDKN1A+CDKN1B+CDKN2A+CHEK1+CHEK2+CIAPIN1+CKB+CKM+CKS2+COL18A1+CSNK2A2+CXCR4+Cellularsenescence+DDB2+DDIT4+DDX5+DFNA5+DKK1+DNAdamage+DUSP2+DUSP4+DUSP5+DYRK2+E2F1+ECT2+EDA2R+EGFR+EIF2AK2+ELAVL1+EPHB4+ERBB2+ESR1+EZH2+FAS+FDXR+FEN1+FGF2+FHL2+FOS+FOXM1+GADD45A+GAPDH+GSTP1+GTSE1+H2AFZ+HDAC1+HIC1+HIF1A+HIPK2+HIPK4+HMMR+HNF4A+HOXA11+HSP90AB1+HSPA4+HTATIP2+ICAM1+ID3+IER3+IFI16+IFITM2+IFNA1+IGF1R+IGFBP1+IGFBP7+IL6+IQCB1+ISG15+KAT2B+KLF4+KRT19+KRT8+LATS2+LTF+MAP4+MAP4K4+MAPK1+MAPK14+MAPK8+MAPK9+MCL1+MCTS1+MDM2+MDM4+MGMT+MMP1+MMP13+MMP2+MSH2+MUC1+MYC+MYCN+NCL+NLRC4+NME1+NOTCH1+NOV+NR2C1+NTN1+PADI4+PARK2+PCBP4+PCNA+PDGFRB+PDRG1+PEG3+PERP+PIDD1+PLAUR+POU4F1+PPM1A+PPM1D+PRC1+PRKCA+PRKD1+PRKDC+PRKG1+PRSS50+PSEN1+PSMD10+PTEN+PTGS2+PTTG1+RAD51+RAF1+RAS+RECQL4+RGCC+RGS16+RPRM+RREB1+RRM2B+S100A2+S100A6+SEMA3B+SERPINB5+SERPINF1+SESN2+SFN+SGK1+SIAH1+SIVA1+SLC2A1+SLC2A4+SLC6A6+SOX4+SP7+TCF7L2+TFDP1+TGFA+TGFB1+THBS1+TIAF1+TIGAR+TLR3+TNFRSF10A+TNFRSF10B+TP53AIP1+TP53I13+TP53INP1+VEGFA+VRK1+WWP1+XAF1+YBX1+ZMAT3, data = data1)
y=summary(res.cox.test)
y
#test 30 genes in table 4

#res.cox.test <- coxph(formula=Surv(data1$survival, data1$status) ~ age+APAF1+AURKA+AURKA+BRCA1+CDC20+CHEK1+CKB+CKS2+DDIT4+E2F1+E2F1+ECT2+EZH2+FEN1+FOXM1+FOXM1+GAPDH+HIF1A+HMMR+HSP90AB1+HSP90AB1+MAPK14+MCTS1+MMP2+MUC1+NCL+NLRC4+PLAUR+PRC1+PRC1+PTTG1+RAS+RECQL4+SFN+SIAH1,data = data1)

# multivariate for 30 genes in Table
res.cox.test <- coxph(formula=Surv(data1$survival, data1$status) ~ CKB+MUC1+FOXM1+E2F1+SFN+CKS2+CHEK1+HSP90AB1+RECQL4+PTTG1+AURKA+PRC1+HMMR+FEN1+DDIT4+MMP2+NLRC4+CDC20+RAS+PLAUR+MCTS1+GAPDH+ECT2+EZH2+APAF1+BRCA1+HIF1A+MAPK14+NCL+SIAH1,data = data1)

y=summary(res.cox.test) 
#y$wald["pvalue"]

#y$logtest
#y$waldtest

#write.table(res.cox.test, file="multi_variate_Cox.txt")                                                                 	
                                                                  
#likelihood ratio test
covariates <- c("age","AATF",	"ABCB1",	"ABCC1", "AIFM2",		"APAF1",	"AR",	"ARID3A",	"ATF3",	"ATM",	"ATR",	"AURKA",	"AXIN1",	"Apoptosis",	"BAK1",	"BAX",	"BBC3",	"BCL2",	"BCL3",	"BCL6",	"BDKRB1",	"BNIP3L",	"BRCA1",	"BTG2",	"CALD1",	"CASP8",	"CCNA2",	"CCNB1",	"CCND1",	"CCNG1",	"CD44",	"CD58",	"CD59",	"CD82",	"CDC20",	"CDC25A",	"CDK2",	"CDK4",	"CDK5",	"CDK9",	"CDKN1A",	"CDKN1B",	"CDKN2A",	"CHEK1",	"CHEK2",	"CIAPIN1",	"CKB", "CKM",	"CKS2",	"COL18A1",	"CSNK2A2",	"CXCR4",	"Cellularsenescence",	"DDB2",	"DDIT4",	"DDX5",	"DFNA5",	"DKK1",	"DNAdamage",	"DUSP2",	"DUSP4",	"DUSP5",	"DYRK2",	"E2F1",	"ECT2",	"EDA2R",	"EGFR",	"EIF2AK2",	"ELAVL1",	"EPHB4",	"ERBB2",	"ESR1",	"EZH2",	"FAS",	"FDXR",	"FEN1",	"FGF2",	"FHL2",	"FOS",	"FOXM1",	"GADD45A",	"GAPDH",	"GSTP1",	"GTSE1",	"H2AFZ",	"HDAC1",	"HIC1",	"HIF1A",	"HIPK2",	"HIPK4",	"HMMR",	"HNF4A",	"HOXA11",	"HSP90AB1",	"HSPA4",	"HTATIP2",	"ICAM1",	"ID3",	"IER3",	"IFI16",	"IFITM2",	"IFNA1",	"IGF1R",	"IGFBP1",	"IGFBP7",	"IL6",	"IQCB1",	"ISG15",	"KAT2B",	"KLF4",	"KRT19",	"KRT8",	"LATS2",	"LTF",	"MAP4",	"MAP4K4",	"MAPK1",	"MAPK14","MAPK8",	"MAPK9",	"MCL1",	"MCTS1",	"MDM2",	"MDM4",	"MGMT",	"MMP1",	"MMP13",	"MMP2",	"MSH2",	"MUC1",	"MYC",	"MYCN",	"NCL",	"NLRC4",	"NME1",	"NOTCH1",	"NOV",	"NR2C1",	"NTN1",	"PADI4",	"PARK2",	"PCBP4",	"PCNA",	"PDGFRB",	"PDRG1",	"PEG3",	"PERP",	"PIDD1",	"PLAUR",	"POU4F1",	"PPM1A",	"PPM1D",	"PRC1",	"PRKCA",	"PRKD1",	"PRKDC",	"PRKG1",	"PRSS50",	"PSEN1",	"PSMD10",	"PTEN",	"PTGS2",	"PTTG1",	"RAD51",	"RAF1",	"RAS",	"RECQL4",	"RGCC",	"RGS16",	"RPRM",	"RREB1",	"RRM2B",	"S100A2",	"S100A6",	"SEMA3B",	"SERPINB5",	"SERPINF1",	"SESN2",	"SFN",	"SGK1",	"SIAH1",	"SIVA1",	"SLC2A1",	"SLC2A4",	"SLC6A6",	"SOX4",	"SP7",	"TCF7L2",	"TFDP1",	"TGFA",	"TGFB1",	"THBS1",	"TIAF1",	"TIGAR",	"TLR3",	"TNFRSF10A",	"TNFRSF10B","TP53AIP1",	"TP53I13",	"TP53INP1",	"VEGFA",	"VRK1",	"WWP1",	"XAF1",	"YBX1",	"ZMAT3")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(data1$survival, data1$status)~', x)))

#univ_formulas <- sapply(covariates,
#                        function(x) as.formula(paste('Surv(data2$survival, data2$status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data1)})
#univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data2)})


# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$logtest["pvalue"], digits=3)
                         logtest<-signif(x$logtest["test"], digits=3)
                         beta<-signif(x$coef[1], digits=3);#coeficient beta
                         HR <-signif(x$coef[2], digits=3);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, logtest, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "Likelihood.ratio.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

#write.csv(as.data.frame(res),file="Univariate_Cox_regression__for_65_samples.csv")
#write.csv(as.data.frame(res),file="Univariate_Cox_regression__for_71_samples_updated_01082018_log2.csv")
write.csv(as.data.frame(res),file="Univariate_Cox_regression__for_71_samples_by_Likelihood_ratio_test_by100.csv")














#wald test result                                                                  
                                                                  

#res.cox.test <- coxph(formula=Surv(data1$survival, data1$status) ~ age+AATF+ABCB1+ABCC1+AIFM2,data=data1)
#res.cox.test
#write.table(res.cox.test$score, file="testcox.txt",quote= FALSE,sep="")

#res.cox.test_1 <- coxph(formula=as.formula(paste('Surv(data1$survival, data1$status)~', x)), data = data1)
#res.cox.test_1
#res.cox.test_1 <- coxph(formula=Surv(data2$survival, data2$status) ~ age + AATF , data = data2)
#write.csv(as.data.frame(res.cox.test_1),file="Multi_variate_Cox_regression.csv")

#covariates <- c("age","stage","AATF",	"ABCB1",	"ABCC1", "AIFM2",		"APAF1",	"AR",	"ARID3A",	"ATF3",	"ATM",	"ATR",	"AURKA",	"AXIN1",	"Apoptosis",	"BAK1",	"BAX",	"BBC3",	"BCL2",	"BCL3",	"BCL6",	"BDKRB1",	"BNIP3L",	"BRCA1",	"BTG2",	"CALD1",	"CASP8",	"CCNA2",	"CCNB1",	"CCND1",	"CCNG1",	"CD44",	"CD58",	"CD59",	"CD82",	"CDC20",	"CDC25A",	"CDK2",	"CDK4",	"CDK5",	"CDK9",	"CDKN1A",	"CDKN1B",	"CDKN2A",	"CHEK1",	"CHEK2",	"CIAPIN1",	"CKB", "CKM",	"CKS2",	"COL18A1",	"CSNK2A2",	"CXCR4",	"Cellularsenescence",	"DDB2",	"DDIT4",	"DDX5",	"DFNA5",	"DKK1",	"DNAdamage",	"DUSP2",	"DUSP4",	"DUSP5",	"DYRK2",	"E2F1",	"ECT2",	"EDA2R",	"EGFR",	"EIF2AK2",	"ELAVL1",	"EPHB4",	"ERBB2",	"ESR1",	"EZH2",	"FAS",	"FDXR",	"FEN1",	"FGF2",	"FHL2",	"FOS",	"FOXM1",	"GADD45A",	"GAPDH",	"GSTP1",	"GTSE1",	"H2AFZ",	"HDAC1",	"HIC1",	"HIF1A",	"HIPK2",	"HIPK4",	"HMMR",	"HNF4A",	"HOXA11",	"HSP90AB1",	"HSPA4",	"HTATIP2",	"ICAM1",	"ID3",	"IER3",	"IFI16",	"IFITM2",	"IFNA1",	"IGF1R",	"IGFBP1",	"IGFBP7",	"IL6",	"IQCB1",	"ISG15",	"KAT2B",	"KLF4",	"KRT19",	"KRT8",	"LATS2",	"LTF",	"MAP4",	"MAP4K4",	"MAPK1",	"MAPK14","MAPK8",	"MAPK9",	"MCL1",	"MCTS1",	"MDM2",	"MDM4",	"MGMT",	"MMP1",	"MMP13",	"MMP2",	"MSH2",	"MUC1",	"MYC",	"MYCN",	"NCL",	"NLRC4",	"NME1",	"NOTCH1",	"NOV",	"NR2C1",	"NTN1",	"PADI4",	"PARK2",	"PCBP4",	"PCNA",	"PDGFRB",	"PDRG1",	"PEG3",	"PERP",	"PIDD1",	"PLAUR",	"POU4F1",	"PPM1A",	"PPM1D",	"PRC1",	"PRKCA",	"PRKD1",	"PRKDC",	"PRKG1",	"PRSS50",	"PSEN1",	"PSMD10",	"PTEN",	"PTGS2",	"PTTG1",	"RAD51",	"RAF1",	"RAS",	"RECQL4",	"RGCC",	"RGS16",	"RPRM",	"RREB1",	"RRM2B",	"S100A2",	"S100A6",	"SEMA3B",	"SERPINB5",	"SERPINF1",	"SESN2",	"SFN",	"SGK1",	"SIAH1",	"SIVA1",	"SLC2A1",	"SLC2A4",	"SLC6A6",	"SOX4",	"SP7",	"TCF7L2",	"TFDP1",	"TGFA",	"TGFB1",	"THBS1",	"TIAF1",	"TIGAR",	"TLR3",	"TNFRSF10A",	"TNFRSF10B","TP53AIP1",	"TP53I13",	"TP53INP1",	"VEGFA",	"VRK1",	"WWP1",	"XAF1",	"YBX1",	"ZMAT3")

covariates <- c("age","AATF",	"ABCB1",	"ABCC1", "AIFM2",		"APAF1",	"AR",	"ARID3A",	"ATF3",	"ATM",	"ATR",	"AURKA",	"AXIN1",	"Apoptosis",	"BAK1",	"BAX",	"BBC3",	"BCL2",	"BCL3",	"BCL6",	"BDKRB1",	"BNIP3L",	"BRCA1",	"BTG2",	"CALD1",	"CASP8",	"CCNA2",	"CCNB1",	"CCND1",	"CCNG1",	"CD44",	"CD58",	"CD59",	"CD82",	"CDC20",	"CDC25A",	"CDK2",	"CDK4",	"CDK5",	"CDK9",	"CDKN1A",	"CDKN1B",	"CDKN2A",	"CHEK1",	"CHEK2",	"CIAPIN1",	"CKB", "CKM",	"CKS2",	"COL18A1",	"CSNK2A2",	"CXCR4",	"Cellularsenescence",	"DDB2",	"DDIT4",	"DDX5",	"DFNA5",	"DKK1",	"DNAdamage",	"DUSP2",	"DUSP4",	"DUSP5",	"DYRK2",	"E2F1",	"ECT2",	"EDA2R",	"EGFR",	"EIF2AK2",	"ELAVL1",	"EPHB4",	"ERBB2",	"ESR1",	"EZH2",	"FAS",	"FDXR",	"FEN1",	"FGF2",	"FHL2",	"FOS",	"FOXM1",	"GADD45A",	"GAPDH",	"GSTP1",	"GTSE1",	"H2AFZ",	"HDAC1",	"HIC1",	"HIF1A",	"HIPK2",	"HIPK4",	"HMMR",	"HNF4A",	"HOXA11",	"HSP90AB1",	"HSPA4",	"HTATIP2",	"ICAM1",	"ID3",	"IER3",	"IFI16",	"IFITM2",	"IFNA1",	"IGF1R",	"IGFBP1",	"IGFBP7",	"IL6",	"IQCB1",	"ISG15",	"KAT2B",	"KLF4",	"KRT19",	"KRT8",	"LATS2",	"LTF",	"MAP4",	"MAP4K4",	"MAPK1",	"MAPK14","MAPK8",	"MAPK9",	"MCL1",	"MCTS1",	"MDM2",	"MDM4",	"MGMT",	"MMP1",	"MMP13",	"MMP2",	"MSH2",	"MUC1",	"MYC",	"MYCN",	"NCL",	"NLRC4",	"NME1",	"NOTCH1",	"NOV",	"NR2C1",	"NTN1",	"PADI4",	"PARK2",	"PCBP4",	"PCNA",	"PDGFRB",	"PDRG1",	"PEG3",	"PERP",	"PIDD1",	"PLAUR",	"POU4F1",	"PPM1A",	"PPM1D",	"PRC1",	"PRKCA",	"PRKD1",	"PRKDC",	"PRKG1",	"PRSS50",	"PSEN1",	"PSMD10",	"PTEN",	"PTGS2",	"PTTG1",	"RAD51",	"RAF1",	"RAS",	"RECQL4",	"RGCC",	"RGS16",	"RPRM",	"RREB1",	"RRM2B",	"S100A2",	"S100A6",	"SEMA3B",	"SERPINB5",	"SERPINF1",	"SESN2",	"SFN",	"SGK1",	"SIAH1",	"SIVA1",	"SLC2A1",	"SLC2A4",	"SLC6A6",	"SOX4",	"SP7",	"TCF7L2",	"TFDP1",	"TGFA",	"TGFB1",	"THBS1",	"TIAF1",	"TIGAR",	"TLR3",	"TNFRSF10A",	"TNFRSF10B","TP53AIP1",	"TP53I13",	"TP53INP1",	"VEGFA",	"VRK1",	"WWP1",	"XAF1",	"YBX1",	"ZMAT3")

univ_formulas <- sapply(covariates,
                       function(x) as.formula(paste('Surv(data1$survival, data1$status)~', x)))

#univ_formulas <- sapply(covariates,
#                        function(x) as.formula(paste('Surv(data2$survival, data2$status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data1)})
#univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data2)})


# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         wald.test<-signif(x$wald["test"], digits=3)
                         beta<-signif(x$coef[1], digits=3);#coeficient beta
                         HR <-signif(x$coef[2], digits=3);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

#write.csv(as.data.frame(res),file="Univariate_Cox_regression__for_65_samples.csv")
#write.csv(as.data.frame(res),file="Univariate_Cox_regression__for_71_samples_updated_01082018_log2.csv")
write.csv(as.data.frame(res),file="Univariate_Cox_regression__for_71_samples_wald_test_by100.csv")


#find gene included in the table 4

input<-read.table(file="associate_gene_in_table_4.txt",header=FALSE, sep="")
gene_1=as.character(unique(input[,1]))
data2=as.data.frame(res)
output=data2[rownames(data2)%in%gene_1,]
#write.csv(output,file="Univariate_Cox_regression__for_gene_in_table_4_updated_01082018_log2.csv")
#write.csv(output,file="Univariate_Cox_regression__for_gene_in_table_4_log_likelyhood_test.csv")
#write.csv(output,file="Univariate_Cox_regression__for_gene_in_table_4_log_likelyhood_test_by100.csv")

write.csv(output,file="Univariate_Cox_regression__for_gene_in_table_4_wald_test_by100.csv")

