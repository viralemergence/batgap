Again, need to download "MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre" and "taxonomy_mamPhy_5911species.csv" from https://github.com/viralemergence/batgap/blob/master/01_data%20processing/code_to_generate_data.csv.R

## libraries
library(dplyr)
library(ggplot2)
library(metafor)
library(ape)
library(caper)
library(viridis)
library(stringr)
library(reshape2)
library(cowplot)
library(ggrepel)
library(car)

##load infection prevalence datasets from https://github.com/viralemergence/batgap/tree/master/02_dissolve%20data
set_infection_prevalence <- read.csv("set_infection_prevalence.csv")
set_infection_prevalence_alphaonly <- read.csv("set_infection_prevalence_alphaonly.csv")
set_infection_prevalence_betaonly <- read.csv("set_infection_prevalence_betaonly.csv")

set_infection_prevalence_nona <- set_infection_prevalence[-which(is.na(set_infection_prevalence$summer)),,]
set_infection_prevalence_alphaonly_nona <- set_infection_prevalence_alphaonly[-which(is.na(set_infection_prevalence_alphaonly$summer)),,]
set_infection_prevalence_betaonly_nona <- set_infection_prevalence_betaonly[-which(is.na(set_infection_prevalence_betaonly$summer)),,]
colnames(set_infection_prevalence_nona)
set_infection_prevalence_nona <- set_infection_prevalence_nona[,c(1:59,68)]
set_infection_prevalence_alphaonly_nona <- set_infection_prevalence_alphaonly_nona[,c(1:59,68)]
set_infection_prevalence_betaonly_nona <- set_infection_prevalence_betaonly_nona[,c(1:59,68)]

#remove variables w/ less than 3 rows (i.e., liver-2, mystacinidae-1, noctilionidae-2): all data 
nrow(set_infection_prevalence_nona)
set_infection_prevalence_nona <- set_infection_prevalence_nona[-which(set_infection_prevalence_nona$tissue_simplified=="liver"),,]
set_infection_prevalence_nona <- set_infection_prevalence_nona[-which(set_infection_prevalence_nona$family=="MYSTACINIDAE"),,]
set_infection_prevalence_nona <- set_infection_prevalence_nona[-which(set_infection_prevalence_nona$family=="NOCTILIONIDAE"),,]
nrow(set_infection_prevalence_nona)

#remove variables w/ less than 3 rows (i.e., mystacinidae-1, noctilionidae-2): alphacov data only
nrow(set_infection_prevalence_alphaonly_nona)
set_infection_prevalence_alphaonly_nona <- set_infection_prevalence_alphaonly_nona[-which(set_infection_prevalence_alphaonly_nona$family=="MYSTACINIDAE"),,]
set_infection_prevalence_alphaonly_nona <- set_infection_prevalence_alphaonly_nona[-which(set_infection_prevalence_alphaonly_nona$family=="NOCTILIONIDAE"),,]
nrow(set_infection_prevalence_alphaonly_nona)

#remove variables w/ less than 3 rows (i.e., liver-2, noctilionidae-2): betacov data only
nrow(set_infection_prevalence_betaonly_nona)
set_infection_prevalence_betaonly_nona <- set_infection_prevalence_betaonly_nona[-which(set_infection_prevalence_betaonly_nona$tissue_simplified=="liver"),,]
set_infection_prevalence_betaonly_nona <- set_infection_prevalence_betaonly_nona[-which(set_infection_prevalence_betaonly_nona$family=="NOCTILIONIDAE"),,]
nrow(set_infection_prevalence_betaonly_nona)

#change order of factor levels
##detection_type
set_infection_prevalence_nona$detection_type <- as.factor(set_infection_prevalence_nona$detection_type)
set_infection_prevalence_nona$detection_type <- factor(set_infection_prevalence_nona$detection_type, levels=c("single","repeat","pooled"))
set_infection_prevalence_alphaonly_nona$detection_type <- as.factor(set_infection_prevalence_alphaonly_nona$detection_type)
set_infection_prevalence_alphaonly_nona$detection_type <- factor(set_infection_prevalence_alphaonly_nona$detection_type, levels=c("single","repeat","pooled"))
set_infection_prevalence_betaonly_nona$detection_type <- as.factor(set_infection_prevalence_betaonly_nona$detection_type)
set_infection_prevalence_betaonly_nona$detection_type <- factor(set_infection_prevalence_betaonly_nona$detection_type, levels=c("single","repeat","pooled"))
##study_type
set_infection_prevalence_nona$study_type <- as.factor(set_infection_prevalence_nona$study_type)
set_infection_prevalence_alphaonly_nona$study_type <- as.factor(set_infection_prevalence_alphaonly_nona$study_type)
set_infection_prevalence_betaonly_nona$study_type <- as.factor(set_infection_prevalence_betaonly_nona$study_type)
##single.multiple.PCR.method
set_infection_prevalence_nona$single.multiple.PCR.method <- as.factor(set_infection_prevalence_nona$single.multiple.PCR.method)
set_infection_prevalence_alphaonly_nona$single.multiple.PCR.method <- as.factor(set_infection_prevalence_alphaonly_nona$single.multiple.PCR.method)
set_infection_prevalence_betaonly_nona$single.multiple.PCR.method <- as.factor(set_infection_prevalence_betaonly_nona$single.multiple.PCR.method)
set_infection_prevalence_nona$single.multiple.PCR.method <- factor(set_infection_prevalence_nona$single.multiple.PCR.method, levels=c("single","multiple"))
set_infection_prevalence_alphaonly_nona$single.multiple.PCR.method <- factor(set_infection_prevalence_alphaonly_nona$single.multiple.PCR.method, levels=c("single","multiple"))
set_infection_prevalence_betaonly_nona$single.multiple.PCR.method <- factor(set_infection_prevalence_betaonly_nona$single.multiple.PCR.method, levels=c("single","multiple"))
##tissue_simplified
set_infection_prevalence_nona$tissue_simplified <- as.factor(set_infection_prevalence_nona$tissue_simplified)
set_infection_prevalence_alphaonly_nona$tissue_simplified <- as.factor(set_infection_prevalence_alphaonly_nona$tissue_simplified)
set_infection_prevalence_betaonly_nona$tissue_simplified <- as.factor(set_infection_prevalence_betaonly_nona$tissue_simplified)
set_infection_prevalence_nona$tissue_simplified <- factor(set_infection_prevalence_nona$tissue_simplified, levels=c("faecal, rectal, or anal","blood or serum","intestine","lung or respiratory","oropharyngeal","skin swab","urinary","pooled swabs/samples","pooled tissue"))
set_infection_prevalence_alphaonly_nona$tissue_simplified <- factor(set_infection_prevalence_alphaonly_nona$tissue_simplified, levels=c("faecal, rectal, or anal","blood or serum","intestine","lung or respiratory","oropharyngeal","skin swab","urinary","pooled swabs/samples","pooled tissue"))
set_infection_prevalence_betaonly_nona$tissue_simplified <- factor(set_infection_prevalence_betaonly_nona$tissue_simplified, levels=c("faecal, rectal, or anal","blood or serum","intestine","lung or respiratory","oropharyngeal","skin swab","urinary","pooled swabs/samples","pooled tissue"))
##gene_targets_simplified
set_infection_prevalence_nona$gene_targets_simplified <- as.factor(set_infection_prevalence_nona$gene_targets_simplified)
set_infection_prevalence_alphaonly_nona$gene_targets_simplified <- as.factor(set_infection_prevalence_alphaonly_nona$gene_targets_simplified)
set_infection_prevalence_betaonly_nona$gene_targets_simplified <- as.factor(set_infection_prevalence_betaonly_nona$gene_targets_simplified)
set_infection_prevalence_nona$gene_targets_simplified <- factor(set_infection_prevalence_nona$gene_targets_simplified, levels=c("RdRp","Other","RdRp_Other"))
set_infection_prevalence_alphaonly_nona$gene_targets_simplified <- factor(set_infection_prevalence_alphaonly_nona$gene_targets_simplified, levels=c("RdRp","Other","RdRp_Other"))
set_infection_prevalence_betaonly_nona$gene_targets_simplified <- factor(set_infection_prevalence_betaonly_nona$gene_targets_simplified, levels=c("RdRp","Other","RdRp_Other"))

## function for I2 for rma.mv
i2=function(model){
  
  ## metafor site code for I2
  W=diag(1/model$vi)
  X=model.matrix(model)
  P=W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2=100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  I2=round(I2,2)
  
  ## summarize by each variance component
  allI2=100 * model$sigma2 / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  allI2=round(allI2,3)
  return(list(I2=I2,allI2=allI2))
}

## trim tree to species in set
tree <- read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")
## load in taxonomy
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa=taxa[taxa$ord=="CHIROPTERA",]
taxa$tip=taxa$Species_Name
taxa$tip <- as.character(taxa$tip)
taxa$fam <- as.character(taxa$fam)

## trim phylo to bats
taxa$tiplabel <- as.character(taxa$tiplabel)
tree <- keep.tip(tree,taxa$tiplabel)

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))
stree1a=keep.tip(tree,as.character(unique(set_infection_prevalence_nona$species)))
stree2b=keep.tip(tree,as.character(unique(set_infection_prevalence_alphaonly_nona$species)))
stree3c=keep.tip(tree,as.character(unique(set_infection_prevalence_betaonly_nona$species)))

## convert tree to correlation matrix
#all genera
cmatrix1a=vcv.phylo(stree1a,cor=T)
#alpha only
cmatrix2b=vcv.phylo(stree2b,cor=T)
#beta only
cmatrix3c=vcv.phylo(stree3c,cor=T)

## make observation and study-level random effect
#all genera
set_infection_prevalence_nona$observation=factor(1:nrow(set_infection_prevalence_nona))
set_infection_prevalence_nona$study=factor(set_infection_prevalence_nona$title)
#alpha only
set_infection_prevalence_alphaonly_nona$observation=factor(1:nrow(set_infection_prevalence_alphaonly_nona))
set_infection_prevalence_alphaonly_nona$study=factor(set_infection_prevalence_alphaonly_nona$title)
#beta only
set_infection_prevalence_betaonly_nona$observation=factor(1:nrow(set_infection_prevalence_betaonly_nona))
set_infection_prevalence_betaonly_nona$study=factor(set_infection_prevalence_betaonly_nona$title)

## pft in escalc for yi and vi 
#all genera
set_infection_prevalence_nona=data.frame(set_infection_prevalence_nona,escalc(xi=set_infection_prevalence_nona$positives,ni=set_infection_prevalence_nona$sample,measure="PFT"))
#alpha only
set_infection_prevalence_alphaonly_nona=data.frame(set_infection_prevalence_alphaonly_nona,escalc(xi=set_infection_prevalence_alphaonly_nona$positives,ni=set_infection_prevalence_alphaonly_nona$sample,measure="PFT"))
#beta only
set_infection_prevalence_betaonly_nona=data.frame(set_infection_prevalence_betaonly_nona,escalc(xi=set_infection_prevalence_betaonly_nona$positives,ni=set_infection_prevalence_betaonly_nona$sample,measure="PFT"))

## back transform
#all genera
set_infection_prevalence_nona$backtrans=transf.ipft(set_infection_prevalence_nona$yi,set_infection_prevalence_nona$sample)
#alpha only
set_infection_prevalence_alphaonly_nona$backtrans=transf.ipft(set_infection_prevalence_alphaonly_nona$yi,set_infection_prevalence_alphaonly_nona$sample)
#beta only
set_infection_prevalence_betaonly_nona$backtrans=transf.ipft(set_infection_prevalence_betaonly_nona$yi,set_infection_prevalence_betaonly_nona$sample)

## species and phylo effect
#all genera
set_infection_prevalence_nona$phylo=set_infection_prevalence_nona$species
set_infection_prevalence_nona$species=set_infection_prevalence_nona$phylo
#alpha only
set_infection_prevalence_alphaonly_nona$phylo=set_infection_prevalence_alphaonly_nona$species
set_infection_prevalence_alphaonly_nona$species=set_infection_prevalence_alphaonly_nona$phylo
#beta only
set_infection_prevalence_betaonly_nona$phylo=set_infection_prevalence_betaonly_nona$species
set_infection_prevalence_betaonly_nona$species=set_infection_prevalence_betaonly_nona$phylo


#run models

model_all4gene_nomethod=rma.mv(yi=yi,V=vi,
                               random=list(~1|study/observation,~1|species,~1|phylo),
                               R=list(phylo=cmatrix1a),
                               method="REML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prevalence_nona,
                               control=list(optimizer="optim", optmethod="BFGS"))

model_all4genealpha_nomethod=rma.mv(yi=yi,V=vi,
                                    random=list(~1|study/observation,~1|species,~1|phylo),
                                    R=list(phylo=cmatrix2b),
                                    method="REML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prevalence_alphaonly_nona,
                                    control=list(optimizer="optim", optmethod="BFGS"))

model_all4genebeta_nomethod=rma.mv(yi=yi,V=vi,
                                   random=list(~1|study/observation,~1|species,~1|phylo),
                                   R=list(phylo=cmatrix3c),
                                   method="REML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prevalence_betaonly_nona,
                                   control=list(optimizer="optim", optmethod="BFGS"))
