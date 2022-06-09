## clean environment & plots
rm(list=ls()) 
graphics.off()

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

##set working directory 
setwd("~/Documents/GitHub/batgap/data")

##load infection prevalence datasets
set_infection_prevalence <- read.csv("set_infection_prevalence.csv")
set_infection_prevalence_alphaonly <- read.csv("set_infection_prevalence_alphaonly.csv")
set_infection_prevalence_betaonly <- read.csv("set_infection_prevalence_betaonly.csv")

##load datasets for binary analyses
set_other <- read.csv("set_other.csv")
set_other_alphaonly <- read.csv("set_other_alphaonly.csv")
set_other_betaonly <- read.csv("set_other_betaonly.csv")

##plot of methods used in set_other/alpha/beta 
count_all <- set_other[which(set_other$method_specific_simplified!="NA"),,] %>%
  group_by(method_specific_simplified) %>%
  dplyr::summarise(count=n())
count_all
count_all$Distribution <- round((count_all$count/sum(count_all$count))*100,digits=2)
sum(3.86,0.18,0.89,0.36,93.5,1.24)

count_all_alpha <- set_other_alphaonly[which(set_other_alphaonly$method_specific_simplified!="NA"),,] %>%
  group_by(method_specific_simplified) %>%
  dplyr::summarise(count=n())
count_all_alpha
count_all_alpha$Distribution <- round((count_all_alpha$count/sum(count_all_alpha$count))*100,digits=2)
sum(1.48,0.11,0.0,0.0,98.1,0.33)

count_all_beta <- set_other_betaonly[which(set_other_betaonly$method_specific_simplified!="NA"),,] %>%
  group_by(method_specific_simplified) %>%
  dplyr::summarise(count=n())
count_all_beta
count_all_beta$Distribution <- round((count_all_beta$count/sum(count_all_beta$count))*100,digits=2)
sum(4.75,0.23,1.13,0.45,92.0,1.47)

count_df <- data.frame(Method=c("ELISA","Indirect Immunofluorescence Test","Indirect-Binding Luminex Assay","Lateral Flow Immunoassay","PCR","Western Blot","ELISA","Indirect Immunofluorescence Test","Indirect-Binding Luminex Assay","Lateral Flow Immunoassay","PCR","Western Blot","ELISA","Indirect Immunofluorescence Test","Indirect-Binding Luminex Assay","Lateral Flow Immunoassay","PCR","Western Blot"), Distribution= c(3.86,0.18,0.89,0.36,93.5,1.24,1.48,0.11,0.0,0.0,98.1,0.33,4.75,0.23,1.13,0.45,92.0,1.47), Data=c("Any Coronavirus","Any Coronavirus","Any Coronavirus","Any Coronavirus","Any Coronavirus","Any Coronavirus","Alphacoronavirus Only","Alphacoronavirus Only","Alphacoronavirus Only","Alphacoronavirus Only","Alphacoronavirus Only","Alphacoronavirus Only","Betacoronavirus Only","Betacoronavirus Only","Betacoronavirus Only","Betacoronavirus Only","Betacoronavirus Only","Betacoronavirus Only"))
count_df$Data <- as.factor(count_df$Data)
count_df$Data <- factor(count_df$Data, levels=c("Any Coronavirus","Alphacoronavirus Only","Betacoronavirus Only"))

ggplot(data = count_df, aes(x = Method, y = Distribution, fill = Method)) + 
  theme_bw() +
  geom_bar(stat = "identity", width=.5, position = "dodge") +
  facet_grid(. ~ Data) + 
  ggtitle("Distribution of Methods Used to Identify Presence of Coronaviruses") + 
  geom_text(aes(label = paste(round(Distribution,digits=2),"%",sep="")), size=3.5,vjust = -0.2) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Distribution") + scale_fill_manual(values=c("#de2160","#9c0946","#d48302","#006400","#07596b","#40013c"))

p <- ggplot(data = count_df, aes(x = Data, y = Distribution, fill = Method)) + 
  theme_bw() +
  geom_bar(stat = "identity", width=.75, position = "stack") +
  ggtitle("Methods Used to Identify Presence of Coronaviruses") +
  theme(plot.title = element_text(hjust = .5,size=8)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) + ylab("Distribution") + scale_fill_manual(values=c("#de2160","#9c0946","#d48302","#006400","#07596b","#40013c")) +
  theme(legend.title = element_text(size=8),legend.text = element_text(size = 8)) + theme(axis.title = element_text(size = 8)) 


png("~/Desktop/plot2.png",width=5,height=8,units="in",res=600)
p
dev.off()

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

##model with no covariates
model_nomods=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix1),
             method="REML",mods=~1,data=set_infection_prev,
             control=list(optimizer="optim", optmethod="BFGS"))

summary(model_nomods)
i2(model_nomods)

model_alpha_nomods=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix2),
             method="REML",mods=~1,data=set_infection_prev_alphaonly,
             control=list(optimizer="optim", optmethod="BFGS"))

summary(model_alpha_nomods)
i2(model_alpha_nomods)

model_beta_nomods=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix3),
             method="REML",mods=~1,data=set_infection_prev_betaonly,
             control=list(optimizer="optim", optmethod="BFGS"))

summary(model_beta_nomods)
i2(model_beta_nomods)

#models for aic comparison


#datasets for following models/aic
set_infection_prevalence_nona <- set_infection_prevalence[-which(is.na(set_infection_prevalence$summer)),,]
set_infection_prevalence_alphaonly_nona <- set_infection_prevalence_alphaonly[-which(is.na(set_infection_prevalence_alphaonly$summer)),,]
set_infection_prevalence_betaonly_nona <- set_infection_prevalence_betaonly[-which(is.na(set_infection_prevalence_betaonly$summer)),,]
colnames(set_infection_prevalence_nona)
set_infection_prevalence_nona <- set_infection_prevalence_nona[,c(1:59,68)]
set_infection_prevalence_alphaonly_nona <- set_infection_prevalence_alphaonly_nona[,c(1:59,68)]
set_infection_prevalence_betaonly_nona <- set_infection_prevalence_betaonly_nona[,c(1:59,68)]

#infection prevalence analyses
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

#covariates = study_type, family, single.multiple.PCR.method, tissue_simplified/sample_simplified, detection_method, euthanasia_simplified, summer/winter/spring/fall, detection_type, gene_targets_simplified


table(set_infection_prev_nona$detection_type, set_infection_prev_nona$gene_targets_simplified)
#two zeros in study_type vs. family
#no zeros in study_type vs single/multiple/pcr, no zeros in study_type vs season, no zeros in study_type vs. euthanasia
#two zeros in study_type vs. detection_type
#five zeros in study_type vs. tissue_simplified
#one zero in study_type vs sample_simplified
#one zero in study_type vs. detection_method
#one zero in study_type vs. gene_targets_simplified
#three zeros in family vs. single multiple PCR
#lots of zeros in family vs. tissue simplified & family vs. sample_simplified
#lots of zeros in family vs. detection_method / detection_type
#three zeros in family vs. euthanasia
#lots of zeros in family vs. gene_targets_simplified
#lots of zeros in family vs. winter, 1 in family vs. fall, 3 in family vs. summer/spring
#four zeros in single/multiple PCR vs tissue 
#one zero in single/multiple PCR vs sample
#one zero in single/multiple PCR vs detection_method, gene_targets
#no zeros in single/multiple PCR vs detection_type, no zeros in single/multiple PCR vs euthanasia, season
#lots of zeros in tissue vs sample, tissue vs detection_method, tissue vs detection_type, tissue vs gene targets; four in tissue vs euthanasia, three in tissue vs season
#lots of zeros in sample vs detection_method, two in sample vs detection_type, one in sample vs euthanasia, one in sample vs season, lots in sample vs gene 
#one zero in detection method vs euthanasia, two in detection method vs gene, two in detection method vs detection type, one in detection method vs season
#none in euthanasia vs detection type, none in euthanasia vs gene, none in euthanasia vs season
#one in season vs gene, none in season vs detection_type, one in detection_type vs gene

#min glob model = single.multiple.PCR, euthanasia, season

model_all_nocol_a=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + euthanasia_simplified + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))
model_all_nocol_b=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + euthanasia_simplified + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

list <- c("study_type","family","single.multiple.PCR.method","tissue_simplified","sample_simplified","detection_method","euthanasia_simplified","season","detection_type","gene_targets_simplified")
DescTools::CombSet(list, m=2) 

###all cov models
model_all1=rma.mv(yi=yi,V=vi,
                      random=list(~1|study/observation,~1|species,~1|phylo),
                      R=list(phylo=cmatrix1a),
                      method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                      control=list(optimizer="optim", optmethod="BFGS"))

model_all2=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all3=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all4=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all5=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all6=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all7=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all8=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all9=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all10=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all11=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all12=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

###alpha cov only models
model_alpha1=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha2=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha3=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~study_type + single.multiple.PCR.method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha4=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha5=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha6=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha7=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha8=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha9=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix2b),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_alpha10=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix2b),
                   method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_alpha11=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix2b),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_alpha12=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix2b),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

###beta cov only models 
model_beta1=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta2=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta3=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~study_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta4=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta5=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta6=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta7=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta8=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta9=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix3c),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_beta10=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta11=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta12=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

###all cov models w gene targets
model_all1gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all2gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all3gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all4gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="REML",mods=~detection_type + detection_method + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

#model_all4genealpha=rma.mv(yi=yi,V=vi,
                      #random=list(~1|study/observation,~1|species,~1|phylo),
                      #R=list(phylo=cmatrix2b),
                      #method="REML",mods=~detection_type + detection_method + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                      #control=list(optimizer="optim", optmethod="BFGS"))

model_all4genebeta=rma.mv(yi=yi,V=vi,
                      random=list(~1|study/observation,~1|species,~1|phylo),
                      R=list(phylo=cmatrix3c),
                      method="REML",mods=~detection_type + detection_method + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                      control=list(optimizer="optim", optmethod="BFGS"))

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

summary(model_all4genebeta_nomethod)


model_all5gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all6gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all7gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all8gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all9gene=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1a),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all10gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix1a),
                   method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_all11gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix1a),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_all12gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix1a),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

###alpha cov only models w gene
model_alpha1gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha2gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha3gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~study_type + single.multiple.PCR.method + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha4gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha5gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha6gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha7gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~detection_type + single.multiple.PCR.method + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha8gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha9gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix2b),
                    method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_alpha10gene=rma.mv(yi=yi,V=vi,
                     random=list(~1|study/observation,~1|species,~1|phylo),
                     R=list(phylo=cmatrix2b),
                     method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                     control=list(optimizer="optim", optmethod="BFGS"))

model_alpha11gene=rma.mv(yi=yi,V=vi,
                     random=list(~1|study/observation,~1|species,~1|phylo),
                     R=list(phylo=cmatrix2b),
                     method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                     control=list(optimizer="optim", optmethod="BFGS"))

model_alpha12gene=rma.mv(yi=yi,V=vi,
                     random=list(~1|study/observation,~1|species,~1|phylo),
                     R=list(phylo=cmatrix2b),
                     method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_alphaonly_nona,
                     control=list(optimizer="optim", optmethod="BFGS"))

###beta cov only models w gene
model_beta1gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta2gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta3gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~study_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta4gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta5gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta6gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta7gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~detection_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta8gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta9gene=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix3c),
                   method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_beta10gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix3c),
                    method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_beta11gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix3c),
                    method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))

model_beta12gene=rma.mv(yi=yi,V=vi,
                    random=list(~1|study/observation,~1|species,~1|phylo),
                    R=list(phylo=cmatrix3c),
                    method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified, data=set_infection_prev_betaonly_nona,
                    control=list(optimizer="optim", optmethod="BFGS"))


## model comparison table
mods_all=list(model_all1,model_all2,model_all3,model_all4,model_all5,model_all6,model_all7,model_all8,model_all9,model_all10,model_all11,model_all12,model_all1gene,model_all2gene,model_all3gene,model_all4gene,model_all6gene,model_all7gene,model_all8gene,model_all9gene,model_all10gene,model_all11gene,model_all12gene)
mods_alpha=list(model_alpha1,model_alpha2,model_alpha3,model_alpha4,model_alpha5,model_alpha6,model_alpha7,model_alpha8,model_alpha9,model_alpha10,model_alpha11,model_alpha12,model_alpha1gene,model_alpha2gene,model_alpha3gene,model_alpha4gene,model_alpha5gene,model_alpha6gene,model_alpha7gene,model_alpha8gene,model_alpha9gene,model_alpha10gene,model_alpha11gene,model_alpha12gene)
mods_beta=list(model_beta1,model_beta2,model_beta3,model_beta4,model_beta5,model_beta6,model_beta7,model_beta8,model_beta9,model_beta10,model_beta11,model_beta12,model_beta1gene,model_beta2gene,model_beta3gene,model_beta4gene,model_alpha5gene,model_beta6gene,model_beta7gene,model_beta8gene,model_beta9gene,model_beta10gene,model_beta11gene,model_beta12gene)

## extract and save
mdata=data.frame(k=sapply(mods_all,function(x) length(coef(x))), AICc=sapply(mods_all,AIC))
mdata=mdata[order(mdata$AIC,decreasing=F),]
mdata$delta=mdata$AIC-mdata$AIC[1]
mdata$wi=MuMIn::Weights(mdata$AIC)

mdata=data.frame(k=sapply(mods_alpha,function(x) length(coef(x))), AICc=sapply(mods_alpha,AIC))
mdata=mdata[order(mdata$AIC,decreasing=F),]
mdata$delta=mdata$AIC-mdata$AIC[1]
mdata$wi=MuMIn::Weights(mdata$AIC)

mdata=data.frame(k=sapply(mods_beta,function(x) length(coef(x))), AICc=sapply(mods_beta,AIC))
mdata=mdata[order(mdata$AIC,decreasing=F),]
mdata$delta=mdata$AIC-mdata$AIC[1]
mdata$wi=MuMIn::Weights(mdata$AIC)
AIC(model_beta10)

#best all = model_all4
mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring
#best alpha = model_alpha10
mods=~detection_type + study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring
#best beta = model_beta6gene
mods=~detection_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring + gene_targets_simplified




#more alpha models


#more beta models 












library(plyr)
#Make models make sense 
#model_all
pred=predict(model_all,transf=transf.ipft.hm,targs=list(ni=set_infection_prevalence$sample))
pred=data.frame(pred,set_infection_prevalence)

pred$id=with(pred,paste(study_type,single.multiple.PCR.method,tissue_simplified,detection_method))
pred=pred[!duplicated(pred$id),]
pred$`Significance of Predicted Mean` <- ifelse(pred$ci.lb!=0.0000000000,"Significant","Not Significant")
                                                 
firstplot <- ggplot(pred,aes(study_type,pred,color=detection_method,alpha=`Significance of Predicted Mean`))+
  geom_errorbar(aes(min=ci.lb,max=ci.ub),width=0,size=0.5,position=position_dodge(width=1))+
  geom_point(size=3,shape=15,position=position_dodge(width=1))+
  facet_grid(study_type~tissue_simplified)+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_point(data=set_infection_prevalence,aes(y=prevalence),alpha=0.5,size=1) + theme(axis.text.x = element_text(angle=90))+ ylab("Alpha of Predicted Mean")+xlab("Study Type")+
  scale_color_manual(values=c("#000000","#009292","#ff6db6","#490092","#24ff24","#006ddb","#920000","#db6d00")) + scale_alpha_manual(values=c(.5,1)) + theme(strip.text.x = element_text(size = 8)) + theme(strip.text.y = element_text(size = 8))+
  guides(colour = guide_legend(order = 1,"Detection Method"), alpha = guide_legend(order = 2))

#model_alpha
pred_alpha=predict(model_alpha,transf=transf.ipft.hm,targs=list(ni=set_infection_prevalence_alphaonly$sample))
pred_alpha=data.frame(pred_alpha,set_infection_prevalence_alphaonly)

pred_alpha$id=with(pred_alpha,paste(study_type,single.multiple.PCR.method,tissue_simplified,detection_method))
pred_alpha=pred_alpha[!duplicated(pred_alpha$id),]
pred_alpha$`Significance of Predicted Mean` <- ifelse(pred_alpha$ci.lb!=0.0000000000,"Significant","Not Significant")

firstplotalpha <- ggplot(pred_alpha,aes(study_type,pred_alpha,color=detection_method,alpha=`Significance of Predicted Mean`))+
  geom_errorbar(aes(min=ci.lb,max=ci.ub),width=0,size=0.5,position=position_dodge(width=1))+
  geom_point(size=3,shape=15,position=position_dodge(width=1))+
  facet_grid(study_type~tissue_simplified)+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_point(data=set_infection_prevalence_alphaonly,aes(y=prevalence),alpha=0.5,size=1) + theme(axis.text.x = element_text(angle=90))+ ylab("Alpha of Predicted Mean")+xlab("Study Type")+
  scale_color_manual(values=c("#000000","#009292","#ff6db6","#490092","#24ff24","#006ddb","#920000","#db6d00")) + scale_alpha_manual(values=c(.5,1)) + theme(strip.text.x = element_text(size = 8)) + theme(strip.text.y = element_text(size = 8))+
  guides(colour = guide_legend(order = 1), alpha = guide_legend(order = 2))

#model_beta
pred_beta=predict(model_beta,transf=transf.ipft.hm,targs=list(ni=set_infection_prevalence_betaonly$sample))
pred_beta=data.frame(pred_beta,set_infection_prevalence_betaonly)

pred_beta$id=with(pred_beta,paste(study_type,single.multiple.PCR.method,tissue_simplified,detection_method))
pred_beta=pred_beta[!duplicated(pred_beta$id),]
pred_beta$`Significance of Predicted Mean` <- ifelse(pred_beta$ci.lb!=0.0000000000,"Significant","Not Significant")

firstplotbeta <- ggplot(pred_beta,aes(study_type,pred_beta,color=detection_method,alpha=`Significance of Predicted Mean`))+
  geom_errorbar(aes(min=ci.lb,max=ci.ub),width=0,size=0.5,position=position_dodge(width=1))+
  geom_point(size=3,shape=15,position=position_dodge(width=1))+
  facet_grid(study_type~tissue_simplified)+ theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_point(data=set_infection_prevalence_betaonly,aes(y=prevalence),alpha=0.5,size=1) + theme(axis.text.x = element_text(angle=90))+ ylab("Alpha of Predicted Mean")+xlab("Study Type")+
  scale_color_manual(values=c("#000000","#009292","#ff6db6","#490092","#24ff24","#006ddb","#920000","#db6d00")) + scale_alpha_manual(values=c(.5,1)) + theme(strip.text.x = element_text(size = 8)) + theme(strip.text.y = element_text(size = 8))+
  guides(colour = guide_legend(order = 1), alpha = guide_legend(order = 2))
