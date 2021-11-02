
#This script is used to performe multiple comparisons on general linear model (GLM) and display in boxplot.

#################################
install.packages('agricolae')
install.packages('lme4')
install.packages('pbkrtest')
install.packages('lmerTest')

#################################
library(data.table)
library(magrittr)
library(stringr)
library(tidyr)
library(lmerTest)
library(MASS)
library(pbkrtest)
library(lme4)
install.packages('agricolae')

################################
library(agricolae)

data=read.table("gene1.HAP.PCA.boxplot.xls",sep="\t",head=T)

#duncan.test
data$Hap=as.factor(data$Hap)
model=glm(Phe ~ Hap+PC1+PC2+PC3,data=data )
summary(model)
ANOVA=aov(model)
re=duncan.test(ANOVA,"Hap")
#plot(re)

##BOXPLOT

#groups
LABELS=data.frame(as.factor(rownames(re$groups)),as.factor(re$groups$groups))
colnames(LABELS)=c("Hap","group")
LABELS$group=LABELS$group[order(as.numeric(LABELS$Hap))]
LABELS$Hap=LABELS$Hap[order(as.numeric(LABELS$Hap))]

jpeg("gene1.HAP.boxplot.jpg",units = "cm",height = 4*5.3,width = 5*5.3,res = 300)

# boxplot
colr=rainbow(nlevels(data$Hap),alpha =0.9)
a <- boxplot(data$Phe ~ data$Hap , ylim=c(0.9*min(data$Phe) , 1.1*max(data$Phe)),pch =NA, col=colr[as.numeric(LABELS$group)] , xlab="Haplotypes",ylab="Leaf length (cm)" , main="gene1")

# add group number and letters
over <- 0.05*max( a$stats[nrow(a$stats),] )
text( c(1:nlevels(data$Hap)) , a$stats[nrow(a$stats),]+over ,  paste("n = ",table(data$Hap),sep="")  )
text( c(1:nlevels(data$Hap)) , a$stats[nrow(a$stats),]+over+4 , LABELS[,2])

dev.off()

################################################################################################
#multiple comparisonsin agricolae package
re1=duncan.test(model,"Hap")
re1$groups
re2 = LSD.test(model,"Hap")
re2$groups
################################################################################################