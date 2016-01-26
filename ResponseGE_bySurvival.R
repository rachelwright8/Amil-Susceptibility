library("DESeq2") # for DEG analysis 
library("arrayQualityMetrics") # for outlier detection
library("Biobase") #for ExpressionSet (required for arrayQualityMetrics)
setwd() # Set your working directory here

#############   #############   #############   #############   #############
#############              Load data                            #############
#############   #############   #############   #############   #############

countdata = read.table("counts.txt",header=TRUE,row.names=1) 
head(countdata) 
length(countdata[,1])
# 44687 genes mapped

names(countdata)
names(countdata) = sub("*.fq.trim.sam.counts","",names(countdata))
names(countdata) = gsub("[.]","",names(countdata)) 

#----------Total counts?
totalCounts = colSums(countdata)
head(totalCounts)
totalCounts
min(totalCounts) #22,116
mean(totalCounts) #360,641.7
max(totalCounts)  #894,289

## Total counts for response:
names(totalCounts)
rtotalCounts = totalCounts[c(1:55,64:124)]
mean(rtotalCounts) #343,773.5
barplot(rtotalCounts)

#-----------Load design table
coldata = read.csv("samples.csv")
head(coldata)

#----get rid of baseline (response only for now)
coldata = coldata[!is.na(coldata$bac),]
coldata

coldata$geno=as.factor(coldata$geno)

#---- combine V. diazotrophicus and V. owensii treatments into a combined bacterial treatment "Y"

test <- function(x) { 
  if(x == "C") y <- "N"
  if(x == "D") y <- "Y"
  if(x == "O") y <- "Y"
  return(y)
}

coldata$bacYN <- sapply(coldata$bac,test)

coldata

#-----------Match order of sample names and conds
coldata$sam
names(countdata)
coldata = coldata[match(names(countdata), coldata$sam),]

#----------Subset for no lesion
coldataN = coldata[(coldata$lesion == "N" & !is.na(coldata$lesion)),]
head(coldataN)

head(countdata)
countdataN = countdata[names(countdata) %in% coldataN$sam]
head(countdataN)

#-----------Check count data and col data
test = cbind(names(countdataN),as.vector(coldataN$sam))
test

#############   #############   #############   #############   #############
#############              Construct data object                #############
#############   #############   #############   #############   #############

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countdataN,
  colData = coldataN,
  design = ~ surv + bacYN + surv:bacYN)

dds = ddsFullCountTable

save(coldata, countdata, dds, file="dds.Rdata")

#############   #############   #############   #############   #############
#############              Call outliers                       #############
#############   #############   #############   #############   #############

rld=rlogTransformation(dds, blind=TRUE) # This takes awhile
e = ExpressionSet(assay(rld), AnnotatedDataFrame(as.data.frame(colData(rld))))
arrayQualityMetrics(e, intgroup=c("surv", "bacYN"), force=T)

##----Failed tests
fail2=c("AD73", "ND262")


#--Remove outliers
countdataOut=countdataN[!(names(countdataN) %in% c(fail2))]
head(countdataOut)

coldataOut=coldataN[(!(coldataN$sam %in% c(fail2))),]
head(coldataOut)

save(coldataOut, countdataOut, rld, e, file="data4DESeq2_outRemoved.Rdata")

#---------Set base mean minimum
means = apply(countdataOut,1,mean)
table(means>3)
# FALSE  TRUE 
# 30054 14633 

means3 = names(means[means>3])
head(means3)
length(means3)
#14633

countDataFilt = countdataOut[row.names(countdataOut) %in% means3,]
head(countDataFilt)

totalCountsFilt = colSums(countDataFilt)
totalCountsFilt

min(totalCountsFilt) #73,991
max(totalCountsFilt) #845,832
mean(totalCountsFilt) #362,569.5

#-------check sample order
test = cbind(names(countDataFilt),as.vector(coldataOut$sam))
test

#############   #############   #############   #############   #############
#############    Reconstruct data object  (outliers removed)    #############
#############   #############   #############   #############   #############

ddsFullCountTableFilt <- DESeqDataSetFromMatrix(
  countData = countDataFilt,
  colData = coldataOut,
  design = ~ surv + bacYN + surv:bacYN)

ddsFilt = ddsFullCountTableFilt

save(coldataOut, countDataFilt, ddsFilt, file="ddsFilt.Rdata")

#############   #############   #############   #############   #############
#############              Explore data                       #############
#############   #############   #############   #############   #############

colMeans(counts(ddsFilt))

apply(counts(ddsFilt), 2, quantile, 0:10/10)

#############   #############   #############   #############   #############
#############          DESeq                                    #############
#############   #############   #############   #############   #############

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFilt)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# -- replacing outliers and refitting for 9 genes
# -- DESeq argument 'minReplicatesForReplace' = 7 
# -- original counts are preserved in counts(dds)
# estimating dispersions
# fitting model and testing

save(countDataFilt, coldataOut, ddsFullCountTable, ddsFilt, deds, e, file="SurvTreat_results.Rdata")

#############   #############   #############   #############   #############
#############          PCA and Heatmaps                         #############
#############   #############   #############   #############   #############

##---------------Sample distance heatmap
rld = rlogTransformation(ddsFilt) # This takes awhile
head(assay(rld))
hist(assay(rld))
colnames(rld) = paste(coldataOut$ab, coldataOut$bac, coldataOut$geno, coldataOut$rep, sep="")
head(assay(rld))

library(pheatmap)
srld=assay(rld)
quartz()
pheatmap(cor(srld),border_color=NA, main="SampleHeatmap")

#---------Principal component analysis 
rld_pca <- function (rld, intgroup = c("surv", "bacYN"), ntop = 14347, colors=NULL, legendpos="topright", main=paste("PCA for top",ntop,sep=" "), textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=0.5))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}

quartz()
rld_pca(rld, intgroup=c("surv", "bacYN"), xlim=c(-35, 55), ylim=c(-30,25))

#############   #############   #############   #############   #############
#############          Extract results                          #############
#############   #############   #############   #############   #############

#--------Interaction
resInt <- results( deds, independentFiltering = F)
resInt
# log2 fold change (MAP): surv.bacYNY 
# Wald test p-value: surv.bacYNY 

table(resInt$padj < 0.05)
# FALSE 
# 14633 
table( is.na(resInt$padj) )
# FALSE 
# 14633  

#--------Bacterial Challenge
resBac <- results( deds, independentFiltering = F, contrast=c("bacYN", "Y", "N"))
resBac
# log2 fold change (MAP): bacYN Y vs N 
# Wald test p-value: bacYN Y vs N 

table(resBac$padj < 0.05)
# FALSE  TRUE 
# 14207   426 
table( is.na(resBac$padj) )
# FALSE   
# 14633     

#--------Survival
resSurv <- results( deds, independentFiltering = F, name="surv")
resSurv
# log2 fold change (MAP): surv
# Wald test p-value: surv

table(resSurv$padj < 0.05)
# FALSE  TRUE 
# 12569  2064 
table( is.na(resSurv$padj) )
# FALSE   
# 14633    

write.csv( as.data.frame(resInt), file="resultsInteraction.csv" )
write.csv( as.data.frame(resBac), file="resulltsBac.csv" )
write.csv( as.data.frame(resSurv), file="resulltsSurv.csv" )

save(coldataOut, countDataFilt, ddsFilt, deds, resInt, resBac, resSurv, rld, file="results.Rdata")

#load("results.Rdata")

#############   #############   #############   #############   #############
#############          Diagnostics                              #############
#############   #############   #############   #############   #############

#####-------------Dispersions plot
quartz()
plotDispEsts(deds, main="Dispersion Plot Response")

####-----------MA plot
quartz()
plotMA(resInt, ylim = c(-1, 1), main="MA Plot Interaction") 

quartz()
plotMA(resBac, ylim = c(-1, 1), main="MA Plot Bacterial Challenge") 

quartz()
plotMA(resSurv, ylim = c(-1, 1), main="MA Plot Survival") 

#############   #############   #############   #############   #############
#############          Write results for heatmaps               #############
#############   #############   #############   #############   #############

###--------------Get pvals
head(resInt)
valsI=cbind(resInt$pvalue, resInt$padj)
head(valsI)
colnames(valsI)=c("pval.i", "padj.i")
length(valsI[,1])
table(complete.cases(valsI))

head(resBac)
valsBac=cbind(resBac$pvalue, resBac$padj)
head(valsBac)
colnames(valsBac)=c("pval.b", "padj.b")
length(valsBac[,1])
table(complete.cases(valsBac))

head(resSurv)
valsSurv=cbind(resSurv$pvalue, resSurv$padj)
head(valsSurv)
colnames(valsSurv)=c("pval.s", "padj.s")
length(valsSurv[,1])
table(complete.cases(valsSurv))

######-------------Make rlogdata and pvals table
rldpvals=cbind(srld,valsI, valsBac, valsSurv)
head(rldpvals)
dim(rldpvals)
#14633    97
table(complete.cases(rldpvals))
# TRUE 
# 14633 

write.csv(rldpvals, "response_RLDandPVALS.csv", quote=F)

#############   #############   #############   #############   #############
#############             Write results for GO                  #############
#############   #############   #############   #############   #############

head(resSurv)
logs=data.frame(cbind("gene"=row.names(resSurv),"logP"=round(-log(resSurv$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resSurv$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 7229 7404 
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Surv_logP.csv",sep=",")

head(resInt)
logs=data.frame(cbind("gene"=row.names(resInt),"logP"=round(-log(resInt$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resInt$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1     1 
# 7861 6772  
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Int_logP.csv",sep=",")

head(resBac)
logs=data.frame(cbind("gene"=row.names(resBac),"logP"=round(-log(resBac$pvalue+1e-10,10),1)))
logs$logP=as.numeric(as.character(logs$logP))
sign=rep(1,nrow(logs))
sign[resBac$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1     1 
# 7432 7201
logs$logP=logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_Bac_logP.csv",sep=",")

#############   #############   #############   #############   #############
#############                       Venn                        #############
#############   #############   #############   #############   #############

library(VennDiagram)

rldpvals=as.data.frame(rldpvals)

surv05=row.names(rldpvals[rldpvals$padj.s<0.05 & !is.na(rldpvals$padj.s),])
bac05=row.names(rldpvals[rldpvals$padj.b<0.05 & !is.na(rldpvals$padj.b),])
int05=row.names(rldpvals[rldpvals$padj.i<0.05 & !is.na(rldpvals$padj.i),])

degs05=union(surv05,bac05)
degs05=union(degs05,int05)

length(degs05)
#2107

candidates=list("surv"=surv05, "bac"=bac05, "int"=int05)
quartz()
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral2", "forestgreen", "royalblue1"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkgreen", "white", "white", "white", "blue4"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkgreen", "blue4"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.03),
  cat.pos = 1
);
grid.draw(prettyvenn)
