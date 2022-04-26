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
unique(set_infection_prevalence$virus_genus)
##set working directory 
setwd("~/Documents/GitHub/batgap/data")

(length(which(set_infection_prevalence[which(set_infection_prevalence$virus_genus=="alphacoronavirus/betacoronavirus coinfection"),,]$positives==0))/length(which(set_infection_prevalence$virus_genus=="alphacoronavirus/betacoronavirus coinfection")))*100
(length(which(set_infection_prevalence[which(set_infection_prevalence$virus_genus=="alphacoronavirus and betacoronavirus and independent bat coronavirus"),,]$positives==0))/length(which(set_infection_prevalence$virus_genus=="alphacoronavirus and betacoronavirus and independent bat coronavirus")))*100

set_infection_prevalence[which(set_infection_prevalence$virus_genus=="alphacoronavirus/alphacoronavirus coinfection"),,]

set_infection_prevalence[which(set_infection_prevalence$title=="ecoepidemiology and complete genome comparison of different strains of severe acute respiratory syndrome-related rhinolophus bat coronavirus in china reveal bats as a reservoir for acute, self-limiting infection that allows recombination events"),,]

##load binary datasets
set_other <- read.csv("set_other.csv")
set_other_alphaonly <- read.csv("set_other_alphaonly.csv")
set_other_betaonly <- read.csv("set_other_betaonly.csv")

##load infection prevalence datasets
set_infection_prevalence <- read.csv("set_infection_prevalence.csv")
set_infection_prevalence_alphaonly <- read.csv("set_infection_prevalence_alphaonly.csv")
set_infection_prevalence_betaonly <- read.csv("set_infection_prevalence_betaonly.csv")

##count_all is going to be positives and negatives while only alpha and only beta are going to be mostly positive?
count_all <- set_other[which(set_other$method_specific_simplified!="NA"),,] %>%
  group_by(method_specific_simplified) %>%
  summarise(count=n())
count_all
sum(3.752,0.172,0.862,0.345,93.618,1.251)

count_all_alpha <- set_other_alphaonly[which(set_other_alphaonly$method_specific_simplified!="NA"),,] %>%
  group_by(method_specific_simplified) %>%
  summarise(count=n())
count_all_alpha
sum(0.657,98.904,0.439)

count_all_beta <- set_other_betaonly[which(set_other_betaonly$method_specific_simplified!="NA"),,] %>%
  group_by(method_specific_simplified) %>%
  summarise(count=n())
count_all_beta
sum(16.349,0.545,5.450,2.180,69.482,5.994)

count_df <- data.frame(Method=c("ELISA","Indirect Immunofluorescence Test","Indirect-Binding Luminex Assay","Lateral Flow Immunoassay","PCR","Western Blot","ELISA","Indirect Immunofluorescence Test","Indirect-Binding Luminex Assay","Lateral Flow Immunoassay","PCR","Western Blot","ELISA","Indirect Immunofluorescence Test","Indirect-Binding Luminex Assay","Lateral Flow Immunoassay","PCR","Western Blot"), Distribution= c(3.752,0.172,0.862,0.345,93.618,1.251,0.657,0,0,0,98.904,0.439,16.349,0.545,5.450,2.180,69.482,5.994), Data=c("All","All","All","All","All","All","Alpha Only","Alpha Only","Alpha Only","Alpha Only","Alpha Only","Alpha Only","Beta Only","Beta Only","Beta Only","Beta Only","Beta Only","Beta Only"))

ggplot(data = count_df, aes(x = Method, y = Distribution, fill = Method)) + 
  theme_bw() +
  geom_bar(stat = "identity", width=.5, position = "dodge") +
  facet_grid(. ~ Data) + 
  ggtitle("Distribution of Methods Used to Identify Coronaviruses") + 
  geom_text(aes(label = paste(round(Distribution,digits=2),"%",sep="")), size=3.5,vjust = -0.2) +
  theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Distribution") + scale_fill_manual(values=c("#de2160","#9c0946","#d48302","#006400","#07596b","#40013c"))

unique(set_infection_prevalence$single.multiple.PCR.method)

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

#trim tree to species in set
stree_all=keep.tip(tree,as.character(unique(set_infection_prevalence$species)))
stree_alpha=keep.tip(tree,as.character(unique(set_infection_prevalence_alphaonly$species)))
stree_beta=keep.tip(tree,as.character(unique(set_infection_prevalence_betaonly$species)))

## convert tree to correlation matrix
cmatrix_all=vcv.phylo(stree_all,cor=T)
cmatrix_alpha=vcv.phylo(stree_alpha,cor=T)
cmatrix_beta=vcv.phylo(stree_beta,cor=T)

set_infection_prevalence$euthanasia_simplified <- set_infection_prevalence$euthanasia
set_infection_prevalence$euthanasia_simplified <- revalue(set_infection_prevalence$euthanasia_simplified, c("Yes, BUT for previous study or rabies surveillance"="Yes","Yes, for this study"="Yes","Yes, for this study; also rabies submissions + found carcasses were processed"="Yes","Yes, but only a subset for species ID, tissue tropism work"="Yes","Yes, but only a subset (no reason provided)"="Yes","No, but some bats were found dead and processed for study"="No","N/A; experimental infection"="No","N/A (samples from previous study used)"="Yes"))
set_infection_prevalence_alphaonly$euthanasia_simplified <- set_infection_prevalence_alphaonly$euthanasia
set_infection_prevalence_alphaonly$euthanasia_simplified <- revalue(set_infection_prevalence_alphaonly$euthanasia_simplified, c("Yes, BUT for previous study or rabies surveillance"="Yes","Yes, for this study"="Yes","Yes, for this study; also rabies submissions + found carcasses were processed"="Yes","N/A; experimental infection"="No"))
set_infection_prevalence_betaonly$euthanasia_simplified <- set_infection_prevalence_betaonly$euthanasia
set_infection_prevalence_betaonly$euthanasia_simplified <- revalue(set_infection_prevalence_betaonly$euthanasia_simplified, c("Yes, BUT for previous study or rabies surveillance"="Yes","Yes, for this study"="Yes","Yes, but only a subset for species ID, tissue tropism work"="Yes","Yes, but only a subset (no reason provided)"="Yes","N/A (samples from previous study used)"="Yes"))

nrow(set_infection_prevalence)
#complete dataset: seasons have 191 NAs (8.9%)
##66.9% = zero prevalence 

#alphacovs only: seasons have 20 NAs (4.6%)
#model can't run because method_specific_simplified only has one value (PCR)
##only 4.6% = zero prevalence 

#betacovs only: seasons have 6 NAs (2.3%)
##only 10.1% = zero prevalence 

unique(set_infection_prevalence$method_specific_simplified)

#models all covariates
model_allcovar=rma.mv(yi=yi,V=vi,
                 random=list(~1|study/observation,~1|species,~1|phylo),
                 R=list(phylo=cmatrix_all),
                 method="REML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prevalence,
                 control=list(optimizer="optim", optmethod="BFGS"))

model_alphaallcovar=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix_alpha),
                   method="REML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prevalence_alphaonly,
                   control=list(optimizer="optim", optmethod="BFGS"))

model_betaallcovar=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix_beta),
                  method="REML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prevalence_betaonly,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all_lm <- lm(yi ~ study_type + single.multiple.PCR.method + tissue_simplified + detection_method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prevalence)
model_alpha_lm <- lm(yi ~ study_type + single.multiple.PCR.method + tissue_simplified + detection_method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prevalence_alphaonly)
model_beta_lm <- lm(yi ~ study_type + single.multiple.PCR.method + tissue_simplified + detection_method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prevalence_betaonly)

car::vif(model_all_lm)

library(car)
regclass::VIF(model_allcovar)
vif(model_betaallcovar)

##model with no covariates
model_nomods=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix_all),
             method="REML",mods=~1,data=set_infection_prevalence,
             control=list(optimizer="optim", optmethod="BFGS"))

model_alpha_nomods=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix_alpha),
             method="REML",mods=~1,data=set_infection_prevalence_alphaonly,
             control=list(optimizer="optim", optmethod="BFGS"))

model_beta_nomods=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix_beta),
             method="REML",mods=~1,data=set_infection_prevalence_betaonly,
             control=list(optimizer="optim", optmethod="BFGS"))


##model with study type as mod
model_with_studytype=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix_all),
                                method="REML",mods=~study_type,data=set_infection_prevalence,
                                control=list(optimizer="optim", optmethod="BFGS"))
model_with_studytype_alpha=rma.mv(yi=yi,V=vi,
                            random=list(~1|study/observation,~1|species,~1|phylo),
                            R=list(phylo=cmatrix_alpha),
                            method="REML",mods=~study_type,data=set_infection_prevalence_alphaonly,
                            control=list(optimizer="optim", optmethod="BFGS"))
model_with_studytype_beta=rma.mv(yi=yi,V=vi,
                            random=list(~1|study/observation,~1|species,~1|phylo),
                            R=list(phylo=cmatrix_beta),
                            method="REML",mods=~study_type,data=set_infection_prevalence_betaonly,
                            control=list(optimizer="optim", optmethod="BFGS"))

##model with method_specific as mod
model_with_methodspecific=rma.mv(yi=yi,V=vi,
                            random=list(~1|study/observation,~1|species,~1|phylo),
                            R=list(phylo=cmatrix_all),
                            method="REML",mods=~method_specific_simplified,data=set_infection_prevalence,
                            control=list(optimizer="optim", optmethod="BFGS"))
model_with_methodspecific_alpha=rma.mv(yi=yi,V=vi,
                                 random=list(~1|study/observation,~1|species,~1|phylo),
                                 R=list(phylo=cmatrix_alpha),
                                 method="REML",mods=~method_specific_simplified,data=set_infection_prevalence_alphaonly,
                                 control=list(optimizer="optim", optmethod="BFGS"))
model_with_methodspecific_beta=rma.mv(yi=yi,V=vi,
                                 random=list(~1|study/observation,~1|species,~1|phylo),
                                 R=list(phylo=cmatrix_beta),
                                 method="REML",mods=~method_specific_simplified,data=set_infection_prevalence_betaonly,
                                 control=list(optimizer="optim", optmethod="BFGS"))
#NA rows omitted from all and alpha

##model with euthanasia as mod
model_with_euthanasia=rma.mv(yi=yi,V=vi,
                                 random=list(~1|study/observation,~1|species,~1|phylo),
                                 R=list(phylo=cmatrix_all),
                                 method="REML",mods=~euthanasia,data=set_infection_prevalence,
                                 control=list(optimizer="optim", optmethod="BFGS"))
model_with_euthanasia_alpha=rma.mv(yi=yi,V=vi,
                             random=list(~1|study/observation,~1|species,~1|phylo),
                             R=list(phylo=cmatrix_alpha),
                             method="REML",mods=~euthanasia,data=set_infection_prevalence_alphaonly,
                             control=list(optimizer="optim", optmethod="BFGS"))
model_with_euthanasia_beta=rma.mv(yi=yi,V=vi,
                             random=list(~1|study/observation,~1|species,~1|phylo),
                             R=list(phylo=cmatrix_beta),
                             method="REML",mods=~euthanasia,data=set_infection_prevalence_betaonly,
                             control=list(optimizer="optim", optmethod="BFGS"))

##model with bat family + sampling season as mods 
model_fam_season=rma.mv(yi=yi,V=vi,
                        random=list(~1|study/observation,~1|species,~1|phylo),
                        R=list(phylo=cmatrix_all),
                        method="REML",mods=~family + summer + winter + fall + spring, data=set_infection_prevalence,
                        control=list(optimizer="optim", optmethod="BFGS"))

model_fam_season_alpha=rma.mv(yi=yi,V=vi,
                        random=list(~1|study/observation,~1|species,~1|phylo),
                        R=list(phylo=cmatrix_alpha),
                        method="REML",mods=~family + summer + winter + fall + spring, data=set_infection_prevalence_alphaonly,
                        control=list(optimizer="optim", optmethod="BFGS"))
model_fam_season_beta=rma.mv(yi=yi,V=vi,
                        random=list(~1|study/observation,~1|species,~1|phylo),
                        R=list(phylo=cmatrix_beta),
                        method="REML",mods=~family + summer + winter + fall + spring, data=set_infection_prevalence_betaonly,
                        control=list(optimizer="optim", optmethod="BFGS"))
#Redundant predictors dropped from the beta model, NA rows removed from all three 
length(which(is.na(set_infection_prevalence$latitude) | is.na(set_infection_prevalence$longitude)))/nrow(set_infection_prevalence)

##model with bat family + sampling season + lat/longs as mods 
model_fam_season_latlong=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix_all),
                                method="REML",mods=~family + summer + winter + fall + spring + latitude + longitude, data=set_infection_prevalence,
                                control=list(optimizer="optim", optmethod="BFGS"))
model_fam_season_latlong_alpha=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix_alpha),
                                method="REML",mods=~family + summer + winter + fall + spring + latitude + longitude, data=set_infection_prevalence_alphaonly,
                                control=list(optimizer="optim", optmethod="BFGS"))
model_fam_season_latlong_beta=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix_beta),
                                method="REML",mods=~family + summer + winter + fall + spring + latitude + longitude, data=set_infection_prevalence_betaonly,
                                control=list(optimizer="optim", optmethod="BFGS"))
#Redundant predictors dropped from the three models, NA rows removed 

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
