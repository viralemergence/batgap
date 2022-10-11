#Code until line 185 is same as code_for_meta-analysis_models in https://github.com/viralemergence/batgap/blob/master/03_analysis/code_for_meta-analysis_models.R

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
set_infection_prevalence_nona <- set_infection_prevalence_nona[,c(1:59,63)]
set_infection_prevalence_alphaonly_nona <- set_infection_prevalence_alphaonly_nona[,c(1:59,63)]
set_infection_prevalence_betaonly_nona <- set_infection_prevalence_betaonly_nona[,c(1:59,63)]

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
which(substr(taxa$tip,1,11)=="Miniopterus")
taxa$fam[which(substr(taxa$tip,1,11)=="Miniopterus")] <- "MINIOPTERIDAE"
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

##Start here if you've already run models (https://github.com/viralemergence/batgap/tree/master/03-generate_and_run_REML_models)
#coefficient plot
global_betas_allcov=data.frame(beta=model_all4gene_nomethod$beta, ci.lb=model_all4gene_nomethod$ci.lb, ci.ub=model_all4gene_nomethod$ci.ub, pval=model_all4gene_nomethod$pval)
global_betas_alphacov=data.frame(beta=model_all4genealpha_nomethod$beta, ci.lb=model_all4genealpha_nomethod$ci.lb, ci.ub=model_all4genealpha_nomethod$ci.ub, pval=model_all4genealpha_nomethod$pval)
global_betas_betacov=data.frame(beta=model_all4genebeta_nomethod$beta, ci.lb=model_all4genebeta_nomethod$ci.lb, ci.ub=model_all4genebeta_nomethod$ci.ub, pval=model_all4genebeta_nomethod$pval)
global_betas_allcov$coefficient <- rownames(global_betas_allcov)
global_betas_alphacov$coefficient <- rownames(global_betas_alphacov)
global_betas_betacov$coefficient <- rownames(global_betas_betacov)
global_betas_allcov$type <- "alphacoronavirus or betacoronavirus"
global_betas_alphacov$type <- "alphacoronavirus only"
global_betas_betacov$type <- "betacoronavirus only"
global_betas <- rbind(global_betas_allcov, global_betas_alphacov, global_betas_betacov)
global_betas$type <- as.factor(global_betas$type)
global_betas$type <- factor(global_betas$type, levels=c("alphacoronavirus or betacoronavirus","alphacoronavirus only","betacoronavirus only"))

unique(global_betas$coefficient)                                                                                      
global_betas$coefficient[which(global_betas$coefficient=="intrcpt")] <- "intercept"
global_betas$coefficient[which(global_betas$coefficient=="detection_typerepeat")] <- "repeat sampling" 
global_betas$coefficient[which(global_betas$coefficient=="detection_typepooled")] <- "pooled sampling"  
global_betas$coefficient[which(global_betas$coefficient=="study_typelongitudinal")] <- "longitudinal study" 
global_betas$coefficient[which(global_betas$coefficient=="single.multiple.PCR.methodmultiple")] <- "multiple PCR runs" 
global_betas$coefficient[which(global_betas$coefficient=="tissue_simplifiedblood or serum")] <- "blood or serum sample" 
global_betas$coefficient[which(global_betas$coefficient=="tissue_simplifiedintestine")] <- "intestinal sample"
global_betas$coefficient[which(global_betas$coefficient=="tissue_simplifiedlung or respiratory")] <- "lung or respiratory sample"  
global_betas$coefficient[which(global_betas$coefficient=="tissue_simplifiedoropharyngeal")] <- "oropharyngeal sample" 
global_betas$coefficient[which(global_betas$coefficient=="tissue_simplifiedpooled swabs/samples")] <- "pooled swabs/samples" 
global_betas$coefficient[which(global_betas$coefficient=="tissue_simplifiedpooled tissue")] <- "pooled tissue" 
global_betas$coefficient[which(global_betas$coefficient=="tissue_simplifiedskin swab")] <- "skin sample"
global_betas$coefficient[which(global_betas$coefficient=="tissue_simplifiedurinary")] <- "urinary sample"
global_betas$coefficient[which(global_betas$coefficient=="euthanasia_simplifiedYes")] <- "euthanasia used" 
global_betas$coefficient[which(global_betas$coefficient=="familyEMBALLONURIDAE")] <- "Emballonuridae"  
global_betas$coefficient[which(global_betas$coefficient=="familyHIPPOSIDERIDAE")] <- "Hipposideridae" 
global_betas$coefficient[which(global_betas$coefficient=="familyMEGADERMATIDAE")] <- "Megadermatidae" 
global_betas$coefficient[which(global_betas$coefficient=="familyMOLOSSIDAE")] <- "Molossidae" 
global_betas$coefficient[which(global_betas$coefficient=="familyMORMOOPIDAE")] <- "Mormoopidae"
global_betas$coefficient[which(global_betas$coefficient=="familyNYCTERIDAE")] <- "Nycteridae"  
global_betas$coefficient[which(global_betas$coefficient=="familyPHYLLOSTOMIDAE")] <- "Phyllostomidae" 
global_betas$coefficient[which(global_betas$coefficient=="familyPTEROPODIDAE")] <- "Pteropodiae" 
global_betas$coefficient[which(global_betas$coefficient=="familyRHINOLOPHIDAE")] <- "Rhinolophidae" 
global_betas$coefficient[which(global_betas$coefficient=="familyRHINOPOMATIDAE")] <- "Rhinopomatidae"
global_betas$coefficient[which(global_betas$coefficient=="familyVESPERTILIONIDAE")] <- "Vespertilionidae"
global_betas$coefficient[which(global_betas$coefficient=="familyMINIOPTERIDAE")] <- "Miniopteridae"
global_betas$coefficient[which(global_betas$coefficient=="summer")] <- "summer" 
global_betas$coefficient[which(global_betas$coefficient=="winter")] <- "winter"  
global_betas$coefficient[which(global_betas$coefficient=="fall")] <- "fall" 
global_betas$coefficient[which(global_betas$coefficient=="spring")] <- "spring" 
global_betas$coefficient[which(global_betas$coefficient=="gene_targets_simplifiedOther")] <- "Not RdRp" 
global_betas$coefficient[which(global_betas$coefficient=="gene_targets_simplifiedRdRp_Other")] <- "RdRp and other target"
global_betas$coefficient <- as.factor(global_betas$coefficient)
global_betas$coefficient <- factor(global_betas$coefficient, levels=rev(c("intercept","repeat sampling","pooled sampling","longitudinal study", "multiple PCR runs", "blood or serum sample","intestinal sample","skin sample","lung or respiratory sample","urinary sample","oropharyngeal sample","pooled swabs/samples","pooled tissue","euthanasia used","Emballonuridae","Hipposideridae","Megadermatidae","Miniopteridae","Molossidae","Mormoopidae","Nycteridae","Phyllostomidae","Pteropodiae","Rhinolophidae","Rhinopomatidae","Vespertilionidae","fall","winter","spring","summer","Not RdRp","RdRp and other target")))

global_betas$variable <- as.factor(c("intercept","sampling method", "sampling method", "study format","PCR type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","euthanasia use","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","sampling season","sampling season","sampling season","sampling season","gene target","gene target","intercept","sampling method", "sampling method", "study format","PCR type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","euthanasia use","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","sampling season","sampling season","sampling season","sampling season","gene target","gene target","intercept","sampling method", "sampling method", "study format","PCR type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","tissue type","euthanasia use","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","bat family","sampling season","sampling season","sampling season","sampling season","gene target","gene target"))
global_betas$variable <- factor(global_betas$variable, levels=c("intercept","sampling method","study format","PCR type","tissue type","euthanasia use","bat family","sampling season","gene target"))

list <- c()
for(i in 1:96){
  if(global_betas$ci.lb[i] < 0 & global_betas$ci.ub[i] < 0){
    new_element <- "no"
    list <- c(list, new_element)
  }
  else if(global_betas$ci.lb[i] > 0 & global_betas$ci.ub[i] > 0){
    new_element <- "no"
    list <- c(list, new_element)
  }
  else{
    new_element <- "yes"
    list <- c(list, new_element)
  }
}
global_betas$cicrosseszero <- list
global_betas$cicrosseszero <- as.factor(global_betas$cicrosseszero)

library(MetBrewer)
library("scales")
colors=met.brewer(name="Cross",n=9,type="continuous")
figure3 <- ggplot(global_betas)+ geom_hline(yintercept=0,linetype="dashed",colour="gray") + geom_errorbar(aes(x=coefficient,ymin=ci.lb,ymax=ci.ub,color=variable,alpha=cicrosseszero),width=0,size=1.5)+scale_alpha_manual(values=c(1, 0.25),guide="none")+geom_point(aes(x=coefficient,y=beta,color=variable),size=3.5)+coord_flip() + facet_grid(.~type) + theme_bw(base_size = 14) + labs(x=NULL) + ylab("meta-analysis model coefficient and 95% confidence interval") + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + scale_color_manual(values=colors)
