# This script is used to performe weighted correlation network analysis (WGCNA) in dynamic transcriptome data.

#################################################################################
install.packages("BiocManager")
BiocManager::install("WGCNA")
#################################################################################

library(WGCNA)
library(flashClust)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Data input

#setwd("E:/WGCNA/R-one-setp")
data=read.table("R_all_deg.rpkm.xls",head=T,row.name=1,sep="\t")
dataF=as.data.frame(t(data))

# Power soft thresholds

powers=c(1:30)
sft=pickSoftThreshold(dataF,powerVector=powers)

pdf(file = "R.softthreshold.pdf", width = 12, height = 9)

par(mfrow = c(1,2))  
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold(power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main=paste("Scale independence"))  
abline(h=0.85, col = "red")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.9,col="red")                                       
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold(power)",ylab="Mean Connectivity", type="n",main=paste("Mean connectivity"))  
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")

dev.off()

#################################################################################

# WGCNA

net = blockwiseModules(dataF,corType="pearson", 
                             networkType="unsigned",power=17,minModuleSize=30, maxBlockSize = 20000,
                             mergeCutHeight=0.2,saveTOMs=TRUE, reassignThreshold = 0,
                             pamRespectsDendro=FALSE,saveTOMFileBase="R",numericLabels = T,verbose=3)

#################################################################################

# Cluster plotting

moduleLabels = net$colors
moduleColors = net$colors #labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

pdf (file="R.Cluster-Dendrogram.pdf")
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#################################################################################

# Output result 

geneClors=data.frame(rownames(data),moduleColors)
write.table(geneClors,file="R.netcolor.xls",sep="\t")
write.table(MEs,file="R.MEs.xls",row.names = rownames(dataF),col.names = T,sep="\t")
save(MEs, moduleLabels, moduleColors,geneClors, geneTree,file = "MEs.RData")

#################################################################################

# ME plots

m1=t(MEs)
colnames(m1)=c("CK","1h","2h","6h","12h","24h")
name=c("CK","1h","2h","6h","12h","24h")

pdf(file="plot-R-MEs.pdf",height = 10,width = 12)
par(mfrow=c(5, 5),mar=c(4,3.5,2,1.5))

for( i in 1:nrow(m1))
{  
	plot(m1[i,],
       type="l",
       ylim=c(-1,1),
       xaxt="n",yaxt="n",
       xlab = "",
       ylab="Module Eigengene (ME)",
       lwd=2,
       col="blue", 
       cex.lab=0.9,
       mgp=c(2.5, 1, 0)
       
  )
  axis(side=1,at=c(1:6),labels = name,cex.axis=1)
  axis(side=2,at=seq(-1,1),cex.axis=1)
  title(sub=rownames(m1)[i],font.sub=2,cex.sub=1.5,mgp=c(2, 1, 0))
}

dev.off()

#################################################################################

# Hug genes

nGenes = ncol(dataF)
nSamples = nrow(dataF)
datKME=signedKME(dataF, MEs)
datKME.P=t(corPvalueStudent(t(datKME),nGenes))
write.table(datKME,file="R.KME.xls",sep="\t",col.names=T)
write.table(datKME.P,file="R.KME-P.xls",sep="\t",col.names=T)

#################################################################################

# Cytoscape 

TOM = TOMsimilarityFromExpr(dataF, power=17) 
probes=names(dataF)
cyt = exportNetworkToCytoscape(TOM,edgeFile="R-sf17.edge.txt",nodeFile="R-sf17.node.txt",weighted = TRUE,nodeNames=probes) 
save(TOM,file="R_sf17.TOM.Rdata")

#################################################################################