library("DESeq2") # for DEG analysis 
library("arrayQualityMetrics") # for outlier detection
library("Biobase") #for ExpressionSet (required for arrayQualityMetrics)
setwd() # Set your working directory here

#############   #############   #############   #############   #############
#############             Load raw count dat                    #############
#############   #############   #############   #############   #############

countdata = read.table("counts.txt",header=TRUE,row.names=1) 

head(countdata) 
length(countdata[,1])
# 44687 isogroups "genes" mapped

names(countdata) # the names are messy, make them better
names(countdata) = sub("*fq.trim.sam.counts","",names(countdata))
names(countdata) = gsub("[.]","",names(countdata)) 

#---------Subset for baseline corals only
countdataB = countdata[c(56:63,125:132)]
head(countdataB)
names(countdataB)

# put them in order.... (27 and 22 are out of place)
countdataB = countdataB[c(1:9,14,11:12,10,13,15:16)]
names(countdataB)
head(countdataB)

#----------Total counts?
totalCounts = colSums(countdataB)
head(totalCounts)
mean(totalCounts) #482936.1
barplot(totalCounts)

min(totalCounts) #192245
max(totalCounts)  #795293

#-----------Load design table. This table has each sample name, reef of origin, treatment, and survival fraction (out of 18 total fragments)
coldata = read.csv("samples.csv")
head(coldata)
# The genotypes are coded as numbers, but they are really factors. 
coldata$geno = as.factor(coldata$geno) 

#----------Subset for baseline
coldataB = coldata[(is.na(coldata$abrasion)),]
coldataB$sam
coldataB

#---------Set base mean minimum = 3
means = apply(countdataB,1,mean)
table(means>3)
#FALSE  TRUE 
#27906 16781 

means3 = names(means[means>3])
head(means3)
length(means3)
#16781

countDataFilt = countdataB[row.names(countdataB) %in% means3,]
head(countDataFilt)

totalCountsFilt = colSums(countDataFilt)
totalCountsFilt

min(totalCountsFilt) #185,392
max(totalCountsFilt) #761,823
mean(totalCountsFilt) #463,096.3

#-------check sample order to make sure samples in counts and conditions folder are aligned
test = cbind(names(countDataFilt),as.vector(coldataB$sam))
test

#############   #############   #############   #############   #############
#############              Construct data object                #############
#############   #############   #############   #############   #############

ddsFullCountTableSurv <- DESeqDataSetFromMatrix(
  countData = countDataFilt,
  colData = coldataB,
  design = ~ reef + surv)

ddsB = ddsFullCountTableSurv

#############   #############   #############   #############   #############
#############              Call outliers                       #############
#############   #############   #############   #############   #############

rld = rlogTransformation(ddsB, blind=T)
e = ExpressionSet(assay(rld), AnnotatedDataFrame(as.data.frame(colData(rld))))
arrayQualityMetrics(e, intgroup=c("surv"), force=T)

#------------No sample outliers

#############   #############   #############   #############   #############
#############          DESeq                                    #############
#############   #############   #############   #############   #############

#-------------DESeq pipeline in one step: makes large DESeqDataSet
dedsB <- DESeq(ddsB)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing


#############   #############   #############   #############   #############
#############          Sample  Heatmaps                         #############
#############   #############   #############   #############   #############
rld = rlogTransformation(dedsB)
head(assay(rld))
colnames(rld) = paste(coldataB$reef, coldataB$geno, coldataB$rep,sep="")
head(assay(rld))

library(pheatmap)
quartz()
pheatmap(cor(assay(rld)),border_color=NA, main="SampleHeatmap")


############         ############         ############         ############                
############                        PCA                        ############ 
############         ############         ############         ############    

rld_pca <- function (rld, intgroup = c("surv", "reef"), ntop = 16000, colors=NULL, legendpos="topright", main=paste("PCA for top",ntop,sep=" "), textcx=1, ...) {
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
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  #legend(legendpos, legend=levels(fac), col=colors, pch=20)
}

quartz()
rld_pca(rld, intgroup=c("surv", "reef"), xlim=c(-45, 50), ylim=c(-50, 40))

#############   #############   #############   #############   #############
#############          Extract results                          #############
#############   #############   #############   #############   #############

resSurv = results(dedsB, name = "surv")
resSurv
# log2 fold change (MAP): surv 
# high survival = 1

resReef = results(dedsB, contrast=c("reef", "W", "L"))
resReef
# log2 fold change (MAP): reef W vs L 
# so positive fold change values indicate a gene that is upregulated at Wilkie

save(coldataB, countDataFilt, rld, ddsB, dedsB, resSurv, resReef, file="baselineDEseq2_bySurv.Rdata")
#load("baselineDEseq2_bySurv.Rdata")

table(resSurv$padj < 0.05)
# FALSE  TRUE 
# 16106   668 

table(resReef$padj < 0.05)
# FALSE  TRUE 
# 15741  1033 

table( is.na(resSurv$padj) )
# FALSE  TRUE 
# 16774     7 

table( is.na(resReef$padj) )
# FALSE  TRUE 
# 16774     7 

resSigSurv <- resSurv[ which(resSurv$padj < 0.05 ), ]
head( resSigSurv[ order( resSigSurv$log2FoldChange ), ] )
resSigReef <- resReef[ which(resReef$padj < 0.05 ), ]
head( resSigReef[ order( resSigReef$log2FoldChange ), ] )

#Number of up-regulated genes with FDR < 0.05
nrow(resSurv[resSurv$padj<0.05 & resSurv$log2FoldChange > 0 & !is.na(resSurv$padj),]) #314
nrow(resReef[resSurv$padj<0.05 & resReef$log2FoldChange > 0 & !is.na(resReef$padj),]) #334

#Number of down-regulated genes with FDR < 0.05
nrow(resSurv[resSurv$padj<0.05 & resSurv$log2FoldChange < 0 & !is.na(resSurv$padj),]) #354
nrow(resReef[resSurv$padj<0.05 & resReef$log2FoldChange < 0 & !is.na(resReef$padj),]) #334

write.csv( as.data.frame(resSurv), file="Surv_results_bm3.csv" )
write.csv( as.data.frame(resReef), file="Reef_results_bm3.csv" )

#####-------------Dispersions plot
quartz()
plotDispEsts(dedsB, main="Dispersion Plot Baseline")

####-----------MA plot
maplot <- function (resSurv, thresh=0.05, thresh2=1, labelsig=TRUE, textcx=1, ...) {
  with(resSurv, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(resSurv, padj<thresh, log2FoldChange>thresh2), points(baseMean, log2FoldChange, col="red", pch=20, cex=0.75))
  if (labelsig) {
    require(calibrate)
    with(subset(resSurv, padj<thresh), textxy(baseMean, log2FoldChange, labs=NA, cex=textcx, col=2))
  }
}
quartz()
maplot(resSurv, main="MA Plot by surv")

maplot <- function (resReef, thresh=0.05, thresh2=1, labelsig=TRUE, textcx=1, ...) {
  with(resReef, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(resReef, padj<thresh, log2FoldChange>thresh2), points(baseMean, log2FoldChange, col="red", pch=20, cex=0.75))
  if (labelsig) {
    require(calibrate)
    with(subset(resReef, padj<thresh), textxy(baseMean, log2FoldChange, labs=NA, cex=textcx, col=2))
  }
}
quartz()
maplot(resReef, main="MA Plot by reef")


#############   #############   #############   #############   #############
#############          Write results for heatmaps               #############
#############   #############   #############   #############   #############

###--------------Get pvals
head(resSurv)
dim(resSurv)
#16781     6
rldd = assay(rld)
head(rldd)
rlddFilt = rldd[row.names(rldd) %in% row.names(resSurv),]
head(rlddFilt)
dim(rlddFilt)
#16781    16

pval.s = resSurv$pvalue
padj.s = resSurv$padj
pval.r = resReef$pvalue
padj.r = resReef$padj
vals = cbind(pval.s,padj.s,pval.r,padj.r)
head(vals)
length(vals[,1])
table(complete.cases(vals))
# FALSE  TRUE 
# 7 16774 

######-------------Make rlogdata and pvals table
head(rlddFilt)
rldpvals = cbind(rlddFilt,vals)
head(rldpvals)
dim(rldpvals)
#16781    18
table(complete.cases(rldpvals))
# TRUE 
# 16781 

write.csv(rldpvals, "baselineSurvReef_BM3_RLDandPVALS.csv", quote=F)


#############   #############   #############   #############   #############
#############             Write results for GO                  #############
#############   #############   #############   #############   #############

head(resSurv)
logs = data.frame(cbind("gene"=row.names(resSurv),"logP" = round(-log(resSurv$pvalue+1e-10,10),1)))
logs$logP = as.numeric(as.character(logs$logP))
sign = rep(1,nrow(logs))
sign[resSurv$log2FoldChange<0] = -1  ##change to correct model
table(sign)
#-1     1 
#8648 8133 
logs$logP = logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_surv_bm3_logP.csv",sep=",")

head(resReef)
logs = data.frame(cbind("gene" = row.names(resReef),"logP" = round(-log(res$pvalue+1e-10,10),1)))
logs$logP = as.numeric(as.character(logs$logP))
sign = rep(1,nrow(logs))
sign[resReef$log2FoldChange<0] = -1  ##change to correct model
table(sign)
#-1     1 
#7879 8902
logs$logP = logs$logP*sign
write.table(logs,quote=F,row.names=F,file="GO_reef_bm3_logP.csv",sep=",")


#############   #############   #############   #############   #############
#############                       Venn                        #############
#############   #############   #############   #############   #############
library(VennDiagram)

rldpvals = as.data.frame(rldpvals)

reef05 = row.names(rldpvals[rldpvals$padj.r<0.05 & !is.na(rldpvals$padj.r),])
surv05 = row.names(rldpvals[rldpvals$padj.s<0.05 & !is.na(rldpvals$padj.s),])

degs05 = union(reef05,surv05)

length(degs05)
#1521

candidates=list("reef"=reef05, "survival"=surv05)
quartz()
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral2", "forestgreen"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkgreen"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)

