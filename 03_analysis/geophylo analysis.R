## bat coronavirus gap analysis
## geographic and taxonomic patterns
## danbeck@ou.edu
## last updated 051223

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(plyr)
library(ggplot2)
library(metafor)
library(ape)
library(caper)
library(viridis)
library(stringr)
library(reshape2)
library(patchwork)

## load datasets for geographic and taxonomic analyses
setwd("~/Desktop/batgap/02_dissolve data")
data_all=read.csv("set_other.csv")
data_alpha=read.csv("set_other_alphaonly.csv")
data_beta=read.csv("set_other_betaonly.csv")

## load in Upham phylogeny
setwd("~/Desktop/batgap/01_data processing")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## load in taxonomy
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa=taxa[taxa$ord=="CHIROPTERA",]
taxa$tip=taxa$Species_Name

## trim phylo to bats
tree=keep.tip(tree,taxa$tiplabel)

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## are all bats in tree
setdiff(data_all$species,tree$tip.label)
setdiff(data_alpha$species,tree$tip.label)
setdiff(data_beta$species,tree$tip.label)

## temporary data file
data=data_all

## clean countries to match world map
data$country=revalue(data$country,
                     c("Malaysian Borneo"="Malaysia",
                       "Lao PDR"="Laos",
                       "United States of America"="USA",
                       "United Kingdom"="UK",
                       "East Timor"="Timor-Leste",
                       "South Korea"="Korea"))

## countries
length(unique(data$country))

## load in global map
library(maptools)
wdata=map_data("world")
wdata=wdata[-which(wdata$region=='Antarctica'),]
wdata$country=wdata$region

## fix Korea
wdata$country=revalue(wdata$country,
                      c("South Korea"="Korea",
                        "North Korea"="Korea"))

## find difference
missing=setdiff(unique(data$country),unique(wdata$country))

## trim from data
set=data[!data$country%in%missing,]
setdiff(set$country,wdata$country)
length(unique(set$country))

## aggregate number of studies and number of bats sampled per country
library(tidyr)
sdata=set %>%
  dplyr::select(studies, species, country, state, site, longitude, latitude, sample_year, start_year, sample) %>%
  dplyr::distinct() %>% 
  dplyr::group_by(country) %>%
  dplyr::summarize(tested = sum(sample),
                   studies = dplyr::n_distinct(studies)) 
adata=data.frame(sdata)
rm(sdata)

## load georegion
setwd("~/Desktop/batgap/03_analysis")
geo=read.csv("georegion.csv")

## standardize
geo$country=geo$name
geo=geo[c("country","region","sub.region")]

## fix country
setdiff(adata$country,geo$country)
geo$country=revalue(geo$country,
                    c("United States of America"="USA",
                      "United Kingdom of Great Britain and Northern Ireland"="UK",
                      "Korea (Democratic People's Republic of)"="Korea",
                      "Korea, Republic of"="Korea",
                      "Congo, Democratic Republic of the"="Democratic Republic of the Congo",
                      "Taiwan, Province of China"="Taiwan",
                      "Trinidad and Tobago"="Trinidad",
                      "Viet Nam"="Vietnam"))
setdiff(adata$country,geo$country)

## merge
gdata=merge(geo,adata,by="country",all.x=T)

## binary sampling
gdata$binstudy=ifelse(is.na(gdata$studies),0,1)

## clean
gdata=gdata[!gdata$region=="",]

## GLM for binary sampling
mod1=glm(binstudy~region,data=gdata,family=binomial)

## within sampled GLMs
gdata2=gdata[gdata$binstudy==1,]
mod2=glm(studies~region,data=gdata2,family=poisson)
mod3=glm(tested~region,data=gdata2,family=poisson)

## range
range(gdata2$studies)
range(gdata2$tested)

## Anova
library(car)
Anova(mod1)
Anova(mod2)
Anova(mod3)

## R2
library(performance)
r2_mcfadden(mod1)
r2_mcfadden(mod2)
r2_mcfadden(mod3)

## visreg
library(visreg)
visreg(mod1,"region",scale="response")
visreg(mod2,"region",scale="response")
visreg(mod3,"region",scale="response")

## posthoc
library(emmeans)
s1=emmeans(mod1,list(pairwise~region),level=0.95,adjust="fdr",type="response")
s2=emmeans(mod2,list(pairwise~region),level=0.95,adjust="fdr",type="response")
s3=emmeans(mod3,list(pairwise~region),level=0.95,adjust="fdr",type="response")

## export
setwd("~/Desktop/batgap/04_outputs")
write.csv(data.frame(s1$`pairwise differences of region`),"Table S3.csv")
write.csv(data.frame(s2$`pairwise differences of region`),"Table S4.csv")
write.csv(data.frame(s3$`pairwise differences of region`),"Table S5.csv")

## merge with wdata
cdata=dplyr::left_join(wdata,adata,by="country",copy=T)

## fix test and studied
cdata$ntested=ifelse(is.na(cdata$tested),0,cdata$tested)
cdata$nstudies=ifelse(is.na(cdata$studies),0,cdata$studies)

## visualize studies
p1=ggplot(cdata,aes(long,lat))+
  geom_polygon(aes(group=group,fill=studies),size=0.1,colour='white')+
  theme_void()+
  scale_fill_viridis_c(na.value="grey70",option="cividis")+
  coord_map("gilbert",xlim=c(-180,180))+
  theme(legend.position="bottom")+
  guides(fill=guide_colorbar(title="(a) studies",
                             barwidth = 15))

## visualize bats
p2=ggplot(cdata,aes(long,lat))+
  geom_polygon(aes(group=group,fill=log1p(tested)),size=0.1,colour='white')+
  theme_void()+
  scale_fill_viridis_c(na.value="grey70",option="cividis")+
  coord_map("gilbert",xlim=c(-180,180))+
  theme(legend.position="bottom")+
  guides(fill=guide_colorbar(title=expression(paste("(b) ",log("samples + 1"))),
                             barwidth = 15))

## combine
setwd("~/Desktop/batgap/04_outputs")
png("Figure 2.png",width=10,height=5,units="in",res=600)
p1+p2
dev.off()

## trim tree to species in set
stree=keep.tip(tree,as.character(unique(data_all$species)))
length(stree$tip.label)

## number of studies per species: all data
sdata=sapply(unique(data_all$species),function(x){
  
  ## tim
  spec=data_all[data_all$species==x,]
  
  ## unique studies
  l=length(unique(spec$studies))
  return(l)
})

## combine
sdata=data.frame(studies=sdata,
                 species=unique(data_all$species))

## save
s_all=sdata ## 343 species

## repeat for alpha
sdata=sapply(unique(data_alpha$species),function(x){
  
  ## tim
  spec=data_alpha[data_alpha$species==x,]
  
  ## unique studies
  l=length(unique(spec$studies))
  return(l)
})

## combine
sdata=data.frame(studies_alpha=sdata,
                 species=unique(data_alpha$species))

## save
s_alpha=sdata

## beta
sdata=sapply(unique(data_beta$species),function(x){
  
  ## tim
  spec=data_beta[data_beta$species==x,]
  
  ## unique studies
  l=length(unique(spec$studies))
  return(l)
})

## combine
sdata=data.frame(studies_beta=sdata,
                 species=unique(data_beta$species))

## save
s_beta=sdata

## save as sdata
sdata=merge(s_all,s_alpha,by="species",all=T)
sdata=merge(sdata,s_beta,by="species",all=T)
rm(s_all,s_alpha,s_beta)

## number of samples per species
samples_all=aggregate(sample~species,data=data_all,sum)
samples_alpha=aggregate(sample~species,data=data_alpha,sum)
samples_beta=aggregate(sample~species,data=data_beta,sum)

## combine
names(samples_alpha)=c("species","sample_alpha")
names(samples_beta)=c("species","sample_beta")

## merge
sdata=merge(sdata,samples_all,by="species",all=T)
sdata=merge(sdata,samples_alpha,by="species",all=T)
sdata=merge(sdata,samples_beta,by="species",all=T)
rm(samples_all,samples_alpha,samples_beta)

## all species
adata=data.frame(species=tree$tip.label)

## merge
sdata=merge(adata,sdata,by="species",all=T)
rm(adata)

## no studies
sdata$studies=ifelse(is.na(sdata$studies),0,sdata$studies)
sdata$studies_alpha=ifelse(is.na(sdata$studies_alpha),0,sdata$studies_alpha)
sdata$studies_beta=ifelse(is.na(sdata$studies_beta),0,sdata$studies_beta)

## no bats
sdata$tested=ifelse(is.na(sdata$sample),0,sdata$sample)
sdata$tested_alpha=ifelse(is.na(sdata$sample_alpha),0,sdata$sample_alpha)
sdata$tested_beta=ifelse(is.na(sdata$sample_beta),0,sdata$sample_beta)

## binary
sdata$binstudy=ifelse(sdata$studies==0,0,1)
sdata$binstudy_alpha=ifelse(sdata$studies_alpha==0,0,1)
sdata$binstudy_beta=ifelse(sdata$studies_beta==0,0,1)

## merge with taxa
taxa=taxa[c("Species_Name","tiplabel","gen","fam","clade","tip","species")]
sdata=merge(sdata,taxa,by="species")

## tally family sampling
sdata2=sdata[sdata$binstudy==1,]
f1=table(sdata2$fam)
f2=table(sdata$fam)

## reorganize
f1=data.frame(f1)
f2=data.frame(f2)
names(f1)=c("fam","sampled")
names(f2)=c("fam","total")
fdata=merge(f2,f1,by="fam",all.x=T)
rm(f1,f2)

## fix NA
fdata$sampled=ifelse(is.na(fdata$sampled),0,fdata$sampled)

## check family coverage
table(fdata$sampled>0)

## check species coverage
sum(fdata$sampled)/sum(fdata$total)

## merge
sdata$tip=sdata$species
cdata=comparative.data(phy=tree,data=sdata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)
cdata$data$label=cdata$data$species
cdata$data$Species=cdata$data$species
cdata$data$tip=cdata$data$species

## taxonomy
cdata$data$taxonomy=with(cdata$data,paste(fam,gen,species,sep='; '))

## just sampled bats
sdata=cdata[cdata$data$binstudy==1,]

## D statistic on sampled/not sampled
set.seed(1)
dstat=phylo.d(data=cdata,binvar=binstudy,permut=1000)
dstat

## number of studies and number of bats
hist(log10(sdata$data$studies))
hist(log10(sdata$data$tested))

## transform
sdata$data$lstudies=log10(sdata$data$studies)
sdata$data$ltested=log10(sdata$data$tested)

## range
range(sdata$data$studies)

## pagel's lambda
pmod1=pgls(lstudies~1,data=sdata,lambda="ML") ## lambda = 0.02
pmod2=pgls(ltested~1,data=sdata,lambda="ML") ## lambda = 0.27

## summarize
summary(pmod1)
summary(pmod2)

## Holm rejection procedure
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## this is the line for Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
library(stringr)
pfsum=function(pf){
  
  ## get formula
  chars=as.character(pf$frmla.phylo)
  chars=strsplit(chars," ")[[1]]
  
  ## response
  resp=chars[1]
  
  ## fix
  #resp=ifelse(resp=='cbind(pos, neg)','prevalence',resp)
  resp=ifelse(str_detect(resp,"cbind"),"prevalence",resp)
  
  ## holm
  hp=HolmProcedure(pf)

  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$tip%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}

## GPF for study binary
library(phylofactor)
set.seed(1)
study_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=binstudy~phylo,
             family=binomial,
             algorithm='phylo',nfactors=2,
             min.group.size = 10)

## summarize
HolmProcedure(study_pf)
study_res=pfsum(study_pf)$results

## number of studies
set.seed(1)
nstudies_pf=gpf(Data=sdata$data,tree=sdata$phy,
              frmla.phylo=studies~phylo,
              family=poisson,
              algorithm='phylo',nfactors=5,
              min.group.size = 10)

## summarize
HolmProcedure(nstudies_pf)
nstudies_res=pfsum(nstudies_pf)$results

## number of samples
set.seed(1)
nsamples_pf=gpf(Data=sdata$data,tree=sdata$phy,
                frmla.phylo=sample~phylo,
                family=poisson,
                algorithm='phylo',nfactors=40,
                min.group.size = 10)

## summarize
HolmProcedure(nsamples_pf)
nsamples_res=pfsum(nsamples_pf)$results

## lower/greater
nsamples_res$check=ifelse(nsamples_res$clade>nsamples_res$other,"more","less")
table(nsamples_res$check)

## save trees
library(treeio)
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")
stree=treeio::full_join(as.treedata(sdata$phy),sdata$data,by="label")

## set x max
plus=1
pplus=plus+0.75

## fix palette
AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
afun=function(x){
  a=AlberColours[1:x]
  return(a)
}

## make low and high
pcols=afun(2)

## function to loop and add clades
cadd=function(gg,pf,pmax){
  
  ## ifelse
  if(HolmProcedure(pf)==0){
    gg=gg
  }else{
    
    ## make result
    result=pfsum(pf)$results
    
    ## ifelse 
    if(nrow(result)>pmax){
      result=result[1:pmax,]
    }else{
      result=result
    }
    
    ## set tree
    for(i in 1:nrow(result)){
      
      ## highlight clade
      gg=gg+
        geom_hilight(node=result$node[i],
                     alpha=0.25,
                     fill=ifelse(result$clade>
                                   result$other,pcols[2],pcols[1])[i])+
        
        ## add label
        geom_cladelabel(node = result$node[i], 
                        label = result$factor[i], 
                        offset = 10, 
                        offset.text = 6,
                        fontsize=2)
      
    }
  }
  return(gg)
}

## state pmax
pmax=24

## make base
library(ggtree)
base=ggtree(dtree,size=0.05)
base2=ggtree(stree,size=0.1,branch.length='none',layout="circular")

## binary tree
gg=cadd(base,study_pf,pmax)

## get tree data
tdata=base$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
library(scales)
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend=rescale(tdata$binstudy,c(max(tdata$x),xmax)),
                species=tdata$species)

## fix gg
plot1=gg+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.25,alpha=0.5)

## fix plot1
plot1=gg+
  geom_tippoint(data=dtree,aes(colour=factor(binstudy),
                               alpha=factor(binstudy)),
                shape=15,size=0.5)+
  guides(colour=F,alpha=F)+
  scale_alpha_manual(values=c(0,1))+
  scale_colour_manual(values=c("white","black"))
plot1=plot1+ggtitle("(a) binary studied")

## number of studies
gg=cadd(base2,nstudies_pf,pmax)

## get tree data
tdata=base2$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
library(scales)
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend=rescale(tdata$studies,c(max(tdata$x),xmax)),
                species=tdata$species)

## fix gg
plot2=gg+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.25,alpha=0.5)
plot2=plot2+ggtitle(expression(paste("(b) ",log[10]("studies"))))

## number of samples
#gg=cadd(base2,nsamples_pf,pmax)

## manual
result=pfsum(nsamples_pf)$results

## trim to first x
x=19
result1=result[1:x,]
result2=result[(x+1):24,]

## set gg
gg=base2

## set tree
for(i in 1:nrow(result1)){
  
  ## highlight clade
  gg=gg+
    geom_hilight(node=result1$node[i],
                 alpha=0.25,
                 fill=ifelse(result1$clade>
                               result1$other,pcols[2],pcols[1])[i])+
    
    ## add label
    geom_cladelabel(node = result1$node[i], 
                    label = result1$factor[i], 
                    offset = 10, 
                    offset.text = 6,
                    fontsize=2)
  
}

## next iteration
for(i in 1:nrow(result2)){
  
  ## highlight clade
  gg=gg+
    geom_hilight(node=result2$node[i],
                 alpha=0.25,
                 fill=ifelse(result2$clade>
                               result2$other,pcols[2],pcols[1])[i])+
    
    ## add label
    geom_cladelabel(node = result2$node[i], 
                    label = result2$factor[i], 
                    offset = 10, 
                    offset.text = 6,
                    hjust=0,  vjust=-1,
                    fontsize=2)
  
}

## get tree data
tdata=base2$data

## tips only
tdata=tdata[which(tdata$isTip==T),]

## set x max
xmax=max(tdata$x)+10

## make data frame
library(scales)
samp=data.frame(x=tdata$x,
                y=tdata$y,
                yend=tdata$y,
                xend=rescale(log10(tdata$tested),c(max(tdata$x),xmax)),
                species=tdata$species)

## fix gg
plot3=gg+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.25,alpha=0.5)
plot3=plot3+ggtitle(expression(paste("(c) ",log[10]("samples"))))

## patchwork and export
library(patchwork)
setwd("~/Desktop/batgap/04_outputs")
png("Figure 3.png",width=6,height=6,units="in",res=600)
plot1|(plot2/plot3)+plot_layout(widths=c(2,1))
dev.off()

## supp tables
setwd("~/Desktop/batgap/04_outputs")
write.csv(nstudies_res,"Table S8.csv")
write.csv(nsamples_res,"Table S9.csv")
