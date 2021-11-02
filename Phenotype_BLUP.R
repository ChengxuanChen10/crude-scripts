

#This script is used to performe best linear unbiased prediction (BLUP) for multi-environments phenotype in mixed linear model .

########################################################
#Import libraries:
library(optparse)
library(data.table)
library(magrittr)
library(stringr)
library(lme4)
########################################################

data=read.table("biomass.cat.xls",sep="\t",head=T)
head(data)

data$Set=data$Set %>% as.factor
data$Block=data$Block %>% as.factor
data$Location=data$Location  %>% as.factor
data$Name2 = data$Name2 %>% as.factor
data$dby=data$dby %>% as.numeric

#mixed model
#re.mod1  = lmer(dby~1+Name2+(1|Location)+(1|Location:Set)+(1|Name2:Location),data=data)
#re.mod2 =lmer(dby~1+Name2+(1|Location)+(1|Location:Block)+(1|Name2:Location),data=data)
re.mod3 =lmer(dby~1+Name2+(1|Location)+(1|Location:Set)+(1|Set:Block)+(1|Name2:Location),data=data)

#fix effect
#fixed_eff=fixef(re.mod1)
#Extract fixed effects:
#y1=fixed_eff[1]+fixed_eff[2:length(fixed_eff)]
#fixed_eff=fixef(re.mod2)
#y2=fixed_eff[1]+fixed_eff[2:length(fixed_eff)]
fixed_eff=fixef(re.mod3)
y3=fixed_eff[1]+fixed_eff[2:length(fixed_eff)]


yo=as.data.frame(y3)

rownames(yo)=gsub("Name2","",rownames(yo))

write.table(yo,file = "biomass.blup.xls" ,sep="\t",quote=F,col.names = F)






