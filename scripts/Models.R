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

ggplot(data = count_df, aes(x = Data, y = Distribution, fill = Method)) + 
  theme_bw() +
  geom_bar(stat = "identity", width=.75, position = "stack") +
  ggtitle("Distribution of Methods Used to Identify Presence of Coronaviruses") +
  theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Distribution") + scale_fill_manual(values=c("#de2160","#9c0946","#d48302","#006400","#07596b","#40013c"))

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

##add even more simplified tissue column
#swab, tissue, faeces, sera, pooled (swab/tissue)
set_infection_prev$sample_simplifed <- set_infection_prev$tissue
set_infection_prev$sample_simplifed <- revalue(set_infection_prev$sample_simplifed, c("rectal swab"="swab","enteric content"="tissue","intestines"="tissue","liver"="tissue","lung"="tissue","faeces/anal swab"="pooled","rectal/oral swab"="swab","faecal swab"="swab","respiratory swab"="swab","faecal swab/faeces"="pooled","alimentary specimen"="tissue","respiratory specimen"="tissue","faeces or urine swab"="pooled",
                                                                                      "urine swab"="swab","faecal and/or throat sample"="pooled","rectal/oral swab,serum"="pooled","oral swab"="swab","faeces or anal swab or intestinal content or oropharyngeal swab"="pooled","faeces or rectal tissue/swab"="pooled","throat swab"="swab","anal swab"="swab","roost feces"="faeces","roost feces,faeces"="faeces","pooled tissue (liver/lung/small intestine/brain/kidney/spleen)"="tissue","urine"="tissue",
                                                                                      "faeces or urine"="pooled","fecal swab"="swab","fecal pellets"="faeces","\"saliva, faeces, and urine samples\""="pooled","anal swabs/faecal samples"="pooled","faecal samples"="faeces","serum"="sera","anal swab,nasopharyngeal swabs"="swab","skin swab"="swab","guano"="faeces","oral swab,alimentary specimen"="pooled","intestines,spleen"="tissue","fecal pellets,anal swab"="pooled","anal swab,oral swab"="swab",
                                                                                      "faecal swab/pellet"="pooled","faeces or rectal swab"="pooled","oral swab,rectal swab,serum"="swab","intestines,lung"="tissue","nasopharyngeal swabs,faecal swab"="swab","pharyngeal/anal swab"="swab","Carcass"="tissue","oral swab,\"tissue (lung,liver,spleen)\""="pooled","nasopharyngeal swabs,anal swab"="swab","fecal pellets,rectal swab,oral swab"="pooled","faeces,oral swab"="pooled","oral swab,rectal swab"="swab",
                                                                                      "rectal swab,oral swab"="swab","respiratory/fecal swab"="swab","oral swab,fecal swab"="swab"))
'%ni%' <- Negate("%in%")
list <- unique(set_infection_prev$tissue)
which(list %ni% set_infection_prev_betaonly$tissue)
set_infection_prev_alphaonly$sample_simplifed <- set_infection_prev_alphaonly$tissue
set_infection_prev_alphaonly$sample_simplifed <- revalue(set_infection_prev_alphaonly$sample_simplifed, c("rectal swab"="swab","intestines"="tissue","liver"="tissue","lung"="tissue","faeces/anal swab"="pooled","rectal/oral swab"="swab","faecal swab"="swab","respiratory swab"="swab","faecal swab/faeces"="pooled","alimentary specimen"="tissue","respiratory specimen"="tissue","faeces or urine swab"="pooled",
                                                                                                          "urine swab"="swab","faecal and/or throat sample"="pooled","rectal/oral swab,serum"="pooled","oral swab"="swab","faeces or anal swab or intestinal content or oropharyngeal swab"="pooled","faeces or rectal tissue/swab"="pooled","throat swab"="swab","anal swab"="swab","roost feces"="faeces","roost feces,faeces"="faeces","pooled tissue (liver/lung/small intestine/brain/kidney/spleen)"="tissue","urine"="tissue",
                                                                                                          "faeces or urine"="pooled","fecal swab"="swab","fecal pellets"="faeces","\"saliva, faeces, and urine samples\""="pooled","anal swabs/faecal samples"="pooled","faecal samples"="faeces","serum"="sera","anal swab,nasopharyngeal swabs"="swab","skin swab"="swab","guano"="faeces","oral swab,alimentary specimen"="pooled","intestines,spleen"="tissue","fecal pellets,anal swab"="pooled","anal swab,oral swab"="swab",
                                                                                                          "faecal swab/pellet"="pooled","faeces or rectal swab"="pooled","oral swab,rectal swab,serum"="swab","intestines,lung"="tissue","nasopharyngeal swabs,faecal swab"="swab","pharyngeal/anal swab"="swab","Carcass"="tissue","oral swab,\"tissue (lung,liver,spleen)\""="pooled","fecal pellets,rectal swab,oral swab"="pooled","faeces,oral swab"="pooled","oral swab,rectal swab"="swab",
                                                                                                          "rectal swab,oral swab"="swab","respiratory/fecal swab"="swab","oral swab,fecal swab"="swab"))

set_infection_prev_betaonly$sample_simplifed <- set_infection_prev_betaonly$tissue
set_infection_prev_betaonly$sample_simplifed <- revalue(set_infection_prev_betaonly$sample_simplifed, c("rectal swab"="swab","enteric content"="tissue","intestines"="tissue","liver"="tissue","lung"="tissue","faeces/anal swab"="pooled","rectal/oral swab"="swab","faecal swab"="swab","respiratory swab"="swab","faecal swab/faeces"="pooled","alimentary specimen"="tissue","respiratory specimen"="tissue","faeces or urine swab"="pooled",
                                                                                                        "faecal and/or throat sample"="pooled","rectal/oral swab,serum"="pooled","oral swab"="swab","faeces or anal swab or intestinal content or oropharyngeal swab"="pooled","faeces or rectal tissue/swab"="pooled","throat swab"="swab","anal swab"="swab","roost feces"="faeces","urine"="tissue",
                                                                                                        "faeces or urine"="pooled","fecal swab"="swab","fecal pellets"="faeces","\"saliva, faeces, and urine samples\""="pooled","anal swabs/faecal samples"="pooled","faecal samples"="faeces","serum"="sera","anal swab,nasopharyngeal swabs"="swab","skin swab"="swab","guano"="faeces","oral swab,alimentary specimen"="pooled","intestines,spleen"="tissue","anal swab,oral swab"="swab",
                                                                                                        "faecal swab/pellet"="pooled","faeces or rectal swab"="pooled","oral swab,rectal swab,serum"="swab","intestines,lung"="tissue","nasopharyngeal swabs,faecal swab"="swab","pharyngeal/anal swab"="swab","Carcass"="tissue","oral swab,\"tissue (lung,liver,spleen)\""="pooled","nasopharyngeal swabs,anal swab"="swab","fecal pellets,rectal swab,oral swab"="pooled","faeces,oral swab"="pooled","oral swab,rectal swab"="swab",
                                                                                                        "rectal swab,oral swab"="swab","respiratory/fecal swab"="swab","oral swab,fecal swab"="swab"))

##add more simplified gene target column 
set_infection_prev$gene_targets_simplified <- set_infection_prev$gene_targets
set_infection_prev$gene_targets_simplified <- revalue(set_infection_prev$gene_targets_simplified, c("ORF 1b"="Other","RdRp: nsP12"="RdRp","ORF SARS-Cov Tor2"="Other","UpE, Orf1a, RdRp, N"="RdRp+Other","ORF1b"="Other","RdRp "="RdRp","E, RdRp"="RdRp+Other","RdRp, Pol"="RdRp+Other","pol"="Other","S1"="Other","RdRP"="RdRp","S"="Other","RdRp-1"="RdRp","RdRp-2"="RdRp"))

set_infection_prev_alphaonly$gene_targets_simplified <- set_infection_prev_alphaonly$gene_targets
set_infection_prev_alphaonly$gene_targets_simplified <- revalue(set_infection_prev_alphaonly$gene_targets_simplified, c("RdRp: nsP12"="RdRp","ORF1b"="Other","RdRp "="RdRp","RdRp, Pol"="RdRp+Other","pol"="Other","S1"="Other","RdRP"="RdRp","S"="Other","E, RdRp"="RdRp+Other"))

set_infection_prev_betaonly$gene_targets_simplified <- set_infection_prev_betaonly$gene_targets
set_infection_prev_betaonly$gene_targets_simplified <- revalue(set_infection_prev_betaonly$gene_targets_simplified, c("ORF 1b"="Other","RdRp: nsP12"="RdRp","ORF SARS-Cov Tor2"="Other","UpE, Orf1a, RdRp, N"="RdRp+Other","ORF1b"="Other","RdRp "="RdRp","E, RdRp"="RdRp+Other","RdRp, Pol"="RdRp+Other","pol"="Other","RdRP"="RdRp","RdRp-1"="RdRp","RdRp-2"="RdRp"))

#datasets for following models/aic
set_infection_prev_nona <- set_infection_prev[-which(is.na(set_infection_prev$summer) | set_infection_prev$gene_targets_simplified=="NA"),,]
set_infection_prev_alphaonly_nona <- set_infection_prev_alphaonly[-which(is.na(set_infection_prev_alphaonly$summer) | set_infection_prev_alphaonly$gene_targets_simplified=="NA"),,]
set_infection_prev_betaonly_nona <- set_infection_prev_betaonly[-which(is.na(set_infection_prev_betaonly$summer) | set_infection_prev_betaonly$gene_targets_simplified=="NA"),,]

nrow(set_infection_prev) - nrow(set_infection_prev_nona)
length(which(is.na(set_infection_prev$summer)))
length(which(set_infection_prev$gene_targets_simplified=="NA"))

colnames(set_infection_prev_nona)
#covariates = study_type, family, single.multiple.PCR.method, tissue_simplified/sample_simplified, detection_method, euthanasia_simplified, summer/winter/spring/fall, detection_type, gene_targets_simplified

###all cov models
model_all1=rma.mv(yi=yi,V=vi,
                      random=list(~1|study/observation,~1|species,~1|phylo),
                      R=list(phylo=cmatrix1),
                      method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev,
                      control=list(optimizer="optim", optmethod="BFGS"))

model_all2=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all3=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~study_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all4=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all5=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all6=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + tissue_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all7=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + euthanasia_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all8=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all9=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~study_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all10=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~detection_type + study_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all11=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

model_all12=rma.mv(yi=yi,V=vi,
                  random=list(~1|study/observation,~1|species,~1|phylo),
                  R=list(phylo=cmatrix1),
                  method="ML",mods=~detection_type + single.multiple.PCR.method + sample_simplified + detection_method + family + summer + winter + fall + spring, data=set_infection_prev,
                  control=list(optimizer="optim", optmethod="BFGS"))

###alpha cov only models


###beta cov only models 


model_alphaallcovar=rma.mv(yi=yi,V=vi,
                           random=list(~1|study/observation,~1|species,~1|phylo),
                           R=list(phylo=cmatrix2),
                           method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_alphaonly,
                           control=list(optimizer="optim", optmethod="BFGS"))

model_betaallcovar=rma.mv(yi=yi,V=vi,
                          random=list(~1|study/observation,~1|species,~1|phylo),
                          R=list(phylo=cmatrix3),
                          method="ML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method + euthanasia_simplified + family + summer + winter + fall + spring, data=set_infection_prev_betaonly,
                          control=list(optimizer="optim", optmethod="BFGS"))


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
