## bat coronavirus gap analysis
## 02_species and country
## danbeck@ou.edu

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

## load data manual
setwd("~/Desktop/batgap/data")
data=read.csv("Prevalence-Grid view.csv")

## ## load in Upham phylogeny
setwd("~/Desktop/BeckerLabOU/phylos")
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
setdiff(data$bat_species,tree$tip)

## multiple species?
data$spec=sapply(strsplit(data$bat_species,","),function(x) length(x))

## fix others
data$species=data$bat_species
data$species=revalue(data$species,
                      c("Pipistrellus javanicas"="Pipistrellus javanicus",
                        "Scotphilus kuhlii"="Scotophilus kuhlii",
                        "Rousettus leschenaulti"="Rousettus leschenaultii",
                        "Myotis ricketti"="Myotis pilosus",
                        "Miniopterus schreibersi"="Miniopterus schreibersii",
                        "Rhinolophus pearsoni"="Rhinolophus pearsonii",
                        "Neoromicia cf. zuluensis"="Neoromicia zuluensis",
                        "Pteropus medius"="Pteropus giganteus",
                        "Rousettous aegyptiacus"="Rousettus aegyptiacus",
                        "Taphozous perforates"="Taphozous perforatus",
                        "Hipposideros caffer ruber"="Hipposideros ruber",
                        "Macronycteris commersoni"="Hipposideros commersoni",
                        "Momopterus acetabulosus"="Mormopterus acetabulosus",
                        "Pteropus seychellensis seychellensis"="Pteropus seychellensis",
                        "myotis bechsteinii"="Myotis bechsteinii",
                        "Pipistrellus nanus"="Neoromicia nana",
                        "Miniopterus africanus"="Miniopterus inflatus",
                        "Mesophylla macconnellii"="Mesophylla macconnelli",
                        "Myotis formosus chofukusei"="Myotis formosus",
                        "Pipistrelus abramus"="Pipistrellus abramus",
                        "Myotis capaccini"="Myotis capaccinii",
                        "Myotis oxygnatus"="Myotis blythii",
                        "Epomps buettikoferi"="Epomops buettikoferi",
                        "Austronomus australis"="Tadarida australis",
                        "Nyctophilus major"="Nyctophilus timoriensis",
                        "Nyctiphilus geoffroyi"="Nyctophilus geoffroyi",
                        "Miniopterus orianae oceanensis"="Miniopterus orianae",
                        "Miniopterus orianae bassanii"="Miniopterus orianae",
                        "Haplonicteris fischeri"="Haplonycteris fischeri",
                        "Dobsonia mollucensis"="Dobsonia moluccensis",
                        "Molossus major"="Molossus molossus",
                        "Rhinolophus affinus"="Rhinolophus affinis",
                        "Megaerops kusnotei"="Megaerops kusnotoi",
                        "Lissonycteris angolensis"="Myonycteris angolensis",
                        "Nycteris cf. gambiensis"="Nycteris gambiensis",
                        "Hipposideros cf. ruber"="Hipposideros ruber",
                        "Hipposideros cf. gigas"="Hipposideros gigas",
                        "Pipistrellus minus"="Pipistrellus tenuis",
                        "Rhinolophus borneenis"="Rhinolophus borneensis",
                        "Miniopterus meghrebensis"="Miniopterus maghrebensis",
                        "Rhinoloplus malayanus"="Rhinolophus malayanus",
                        "Craseonycteris thonlongyal"="Craseonycteris thonglongyai",
                        "Coelops frithi"="Coelops frithii",
                        "Rhinolophus malyanus"="Rhinolophus malayanus",
                        "Peking Myotis,Myotis pequinius"="Myotis pequinius",
                        "Miniopterus filiginosus"="Miniopterus fuliginosus",
                        "Hypsugo pulveratus"="Pipistrellus pulveratus",
                        "Tyloncyteris pachypus"="Tylonycteris pachypus",
                        "Rousettus lechenaulti"="Rousettus leschenaultii",
                        "Rhinolophus affiinus"="Rhinolophus affinis",
                        "Vespertilio superans"="Vespertilio sinensis",
                        "Rhinolophus rouxi"="Rhinolophus rouxii",
                        "Myotis formosus flavus"="Myotis formosus",
                        "Pipistrellus taiwanesis"="Pipistrellus taiwanensis",
                        "Myotis fimbriatus taiwanensis"="Myotis fimbriatus",
                        "Barbastella darjelingesis"="Barbastella beijingensis",
                        "Hipposideros armiger terasensis"="Hipposideros armiger",
                        "Coelops frithii formosanus"="Coelops frithii",
                        "Rhinolophus darlingi damarensis"="Rhinolophus darlingi",
                        "Nycticeinops schlieffenii"="Nycticeinops schlieffeni",
                        "Neuromicia helios"="Neoromicia helios",
                        "Neuromicia nana"="Neoromicia nana",
                        "Scotorepens rueppellii"="Scoteanax rueppellii",
                        "Hipposideros terasensis"="Hipposideros armiger",
                        "Stenonycteris lanosus"="Rousettus lanosus",
                        "Myonicteris angolensis"="Myonycteris angolensis",
                        "Mormoops megalohyla"="Mormoops megalophylla",
                        "Artibeus watsoni"="Dermanura watsoni",
                        "Artibeus phaeotis"="Dermanura phaeotis",
                        "Nycteris argae"="Nycteris arge",
                        "Lyssonycteris angolensis"="Myonycteris angolensis",
                        "Pipistrellus capensis"="Neoromicia capensis",
                        "Myotis bocagei"="Myotis bocagii",
                        "Chaerephon pumila"="Chaerephon pumilus",
                        "Hypsugo savii"="Pipistrellus savii",
                        "Myotis mistacinus"="Myotis mystacinus",
                        "Myotis oxygnathus"="Myotis blythii",
                        "Rhinolophus shamelli"="Rhinolophus shameli",
                        "Rhinolophus blythi"="Rhinolophus lepidus",
                        "Rhinolophus monoceros"="Rhinolophus pusillus",
                        "Rhinolophus lobatus"="landeri",
                        "Rhinolophus rhodesiae"="Rhinolophus simulator",
                        "Pteropus seychellensis comorensis"="Pteropus seychellensis",
                        "Chaerephon pusillus"="Chaerephon pumilus",
                        "Taphosius mauritianus"="Taphozous mauritianus",
                        "Chaerephon leucogaster"="Chaerephon pumilus",
                        "Myotis cillolabrum"="Myotis ciliolabrum",
                        "Vampyressa pussilla"="Vampyressa pusilla",
                        "Artibeus literatus"="Artibeus lituratus",
                        "Miniopterus orianae"="Miniopterus oceanensis",
                        "Vampyriscus nymphaea"="Vampyressa nymphaea",
                        "Glossophaga comissarisi"="Glossophaga commissarisi",
                        "Eptesicus furnalis"="Eptesicus furinalis"))

## check
setdiff(data$species,tree$tip)
length(setdiff(data$species,tree$tip))
length(unique(data$species))
## 4 missing from phylogeny
## Myotis rufoniger, Pipistrellus montanus, Pipistrellus taiwanensis, Rhinolophus cornutus

## trim dataset to only phylogeny
set=data[data$species%in%tree$tip.label,]
set$species=factor(set$species)
set$studies=factor(set$Field.25)

## trim tree to species in set
stree=keep.tip(tree,as.character(unique(set$species)))

## convert tree to correlation matrix
cmatrix=vcv.phylo(stree,cor=T)

## make observation and study-level random effect
set$observation=factor(1:nrow(set))
set$study=factor(set$Field.25)

## pft in escalc for yi and vi 
set=data.frame(set,escalc(xi=set$positives,ni=set$sample,measure="PFT"))

## back transform
set$backtrans=transf.ipft(set$yi,set$sample)

## check
plot(set$prevalence,set$backtrans)
abline(0,1)

## species and phylo effect
set$phylo=set$species
set$species=set$phylo

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

## model with no covariates
model=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix),
             method="REML",mods=~1,data=set,
             control=list(optimizer="optim", optmethod="BFGS"))

## aggregate to sampled bats
library(tidyr)
sdata=set %>%
  #filter(!(studies == "highly diversified coronaviruses in neotropical bats")) %>%
   dplyr::select(studies, species, country, state, site, longitude, latitude, sample_year, start_year, sample) %>%
  dplyr::distinct() %>% 
  dplyr::group_by(species) %>%
   dplyr::summarize(tested = sum(sample),
            nstudies = dplyr::n_distinct(studies)) 
adata=data.frame(sdata)

## number of studies
sdata=sapply(unique(set$species),function(x){
  
  ## tim
  spec=set[set$species==x,]
  
  ## unique studies
  l=length(unique(spec$studies))
  return(l)
})

## combine
sdata=data.frame(studies=sdata,
                 species=unique(set$species))

## test merge
test=merge(adata,sdata,by="species")
plot(test$studies,test$nstudies)
cor(test$studies,test$nstudies)
abline(0,1)

## all species
adata=data.frame(species=tree$tip.label)

## merge
sdata=merge(adata,test,by="species",all=T)
rm(adata,test)

## no studies
sdata$studies=ifelse(is.na(sdata$studies),0,sdata$studies)

## no bats
sdata$tested=ifelse(is.na(sdata$tested),0,sdata$tested)

## binary
sdata$binstudy=ifelse(sdata$studies==0,0,1)

## merge with taxa
taxa=taxa[c("Species_Name","tiplabel","gen","fam","clade","tip","species")]
sdata=merge(sdata,taxa,by="species")

## which is weird
sdata$odd=ifelse(tree$tip.label==sdata$species,1,0)

## get betacov ensemble propranks from 2020
setwd("~/Desktop/Fresnel_Jun/Cleaned Files_2020")
batin=read.csv("BatModels_IS.csv")
batout=read.csv("BatModels_OS.csv")

## remove odd rows
batout$fix=as.numeric(batout$Sp)
batout=batout[is.na(batout$fix),]
batout=batout[!is.na(batout$Sp),]

## combine
batold=rbind.data.frame(batin[intersect(names(batin),names(batout))],
                        batout[intersect(names(batin),names(batout))])
batold=batold[c("Sp","Betacov","PropRank","InSample")]

## clean
rm(batin,batout)

## fix names
batold$species=batold$Sp
batold=batold[c("species","PropRank","Betacov")]

## initial merge
test=merge(sdata,batold,by="species",all=T)

## plot proprank and studies
ggplot(test,aes(PropRank,studies))+
  geom_point(alpha=0.5,
             aes(colour=factor(Betacov)))+
  geom_smooth(method="glm",
              method.args=list(family=poisson))+
  theme_bw()+
  guides(colour=F)+
  scale_colour_manual(values=c("grey","red"))+
  labs(x="proportional ranks for bat betacoronavirus hosts",
       y="number of CoV studies per bat species")

## merge
sdata$tip=sdata$species
cdata=comparative.data(phy=tree,data=sdata,names.col=tip,vcv=T,na.omit=F,warn.dropped=T)
cdata$data$label=cdata$data$species
cdata$data$Species=cdata$data$species
cdata$data$tip=cdata$data$species

## taxonomy
cdata$data$taxonomy=with(cdata$data,paste(fam,gen,species,sep='; '))

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
  chars=as.character(pf$frmla.phylo)[-1]
  
  ## response
  resp=chars[1]
  
  ## fix
  #resp=ifelse(resp=='cbind(pos, neg)','prevalence',resp)
  resp=ifelse(str_detect(resp,"cbind"),"prevalence",resp)
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## save model
  model=chars[2]
  
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

## negbin model fcn
model.fcn2 <- function(formula,data,...){
  fit <- tryCatch(MASS::glm.nb(formula,data,...),
                  error=function(e) NA)
  #fit <- do.call
  return(fit)
}

## negbin objective function
obj.fcn2 <- function(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...){
  #if (!'negbin' %in% class(fit) & !'glm' %in% class(fit) & !'lm' %in% class(fit))
  if (!'negbin' %in% class(fit))
  {
    return(0)
  }
  else 
  {
    #fit2 <- MASS::glm.nb(Z.poisson~1,data = fit$model)
    fit$null.deviance-fit$deviance %>% return()
    #fit$twologlik %>% return()
  }
}

## studies: poisson or nb
mean(cdata$data$studies)
var(cdata$data$studies)

## tested: poisson or nb
mean(cdata$data$tested)
var(cdata$data$tested)

## GPF for studies
library(phylofactor)
library(MASS)
set.seed(1)
study_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=studies~phylo,
             model.fcn = model.fcn2,objective.fcn = obj.fcn2,
             cluster.depends='library(MASS)',
             algorithm='phylo',nfactors=15)

## summarize
HolmProcedure(study_pf)
study_res=pfsum(study_pf)$results

## number of bats tested
set.seed(1)
tested_pf=gpf(Data=cdata$data,tree=cdata$phy,
             frmla.phylo=tested~phylo,
             model.fcn = model.fcn2,objective.fcn = obj.fcn2,
             cluster.depends='library(MASS)',
             algorithm='phylo',nfactors=15)

## summarize
HolmProcedure(tested_pf)
tested_res=pfsum(tested_pf)$results

## save tree
library(treeio)
dtree=treeio::full_join(as.treedata(cdata$phy),cdata$data,by="label")

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
cadd=function(gg,pf){
  
  ## ifelse
  if(HolmProcedure(pf)==0){
    gg=gg
  }else{
    
    ## make result
    result=pfsum(pf)$results
    
    ## set tree
    for(i in 1:nrow(result)){
      
      gg=gg+
        geom_hilight(node=result$node[i],
                     alpha=0.25,
                     fill=ifelse(result$clade>
                                   result$other,pcols[2],pcols[1])[i])
    }
  }
  return(gg)
}

## make base
library(ggtree)
base=ggtree(dtree,size=0.1,branch.length='none',layout="circular")
base

## study tree
gg=cadd(base,study_pf)

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
                xend=rescale(tdata$studies,c(max(tdata$x),xmax)),
                species=tdata$species)

## fix gg
gg+
  geom_segment(data=samp,aes(x=x,y=y,xend=xend,yend=yend),size=0.75,alpha=0.5)

## get undersampled species
usamp1=cladeget(study_pf,1)
usamp2=cladeget(study_pf,2)
usamps=c(usamp1,usamp2)

## genera
usamp1=unique(sapply(strsplit(usamp1," "),function(x) x[1]))
usamp2=unique(sapply(strsplit(usamp2," "),function(x) x[1]))

## load in betacov predictions
setwd("~/Desktop/Fresnel_Jun")
new=read.csv("BinaryPredictions.csv")

## ensemble T and betacov 0
new$suspect=ifelse(new$Ensemble==TRUE & new$Betacov==0,1,0)

## just suspect
new=new[new$suspect==1,]

## which are in poorly studies clades
usamps[usamps%in%new$Sp]

## clean countries to match world map
data$country=revalue(data$country,
                  c("Malaysian Borneo"="Malaysia",
                    "Lao PDR"="Laos",
                    "United States of America"="USA",
                    "United Kingdom"="UK",
                    "East Timor"="Timor-Leste",
                    "South Korea"="Korea"))

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

## aggregate number of studies and number of bats sampled per country
library(tidyr)
set$studies=set$Field.25
sdata=set %>%
  #filter(!(studies == "highly diversified coronaviruses in neotropical bats")) %>%
  dplyr::select(studies, species, country, state, site, longitude, latitude, sample_year, start_year, sample) %>%
  distinct() %>% 
  group_by(country) %>%
  dplyr::summarize(tested = sum(sample),
                   studies = n_distinct(studies)) 
adata=data.frame(sdata)
rm(sdata)

## merge with wdata
cdata=left_join(wdata,adata,by="country",copy=T)

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
  guides(fill=guide_colorbar(title="number of studies",
                             barwidth = 15))

## visualize bats
p2=ggplot(cdata,aes(long,lat))+
  geom_polygon(aes(group=group,fill=log10(tested)),size=0.1,colour='white')+
  theme_void()+
  scale_fill_viridis_c(na.value="grey70",option="cividis")+
  coord_map("gilbert",xlim=c(-180,180))+
  theme(legend.position="bottom")+
  guides(fill=guide_colorbar(title="log10 number of bats",
                             barwidth = 15))

## combine
library(patchwork)
p1+p2