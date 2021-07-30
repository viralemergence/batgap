##bat coronavirus gap analysis
##species and country

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

##load data
setwd("~/Desktop")
data <- read.csv("batdata.csv")
##get rid of row #98 (not sure where that number came from upon rereading paper)
data <- data[-c(98),,]

##load in Upham phylogeny 
tree <- read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")

## load in taxonomy
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa=taxa[taxa$ord=="CHIROPTERA",]
taxa$tip=taxa$Species_Name
taxa$tip <- as.character(taxa$tip)

## trim phylo to bats
taxa$tiplabel <- as.character(taxa$tiplabel)
tree <- keep.tip(tree,taxa$tiplabel)

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## are all bats in tree
setdiff(data$bat_species,tree$tip)

## multiple species?
data$bat_species <- as.character(data$bat_species)
data$spec=sapply(strsplit(data$bat_species,","),function(x) length(x))

## fix others
data$species=data$bat_species
data$species=revalue(data$species,
                     c("Eptesicus"="",
                       "Artibeus literatus"="",
                       "Vampyressa pussilla"="",
                       "Myotis cillolabrum"="",
                       "Macronycteris commersoni"="",
                       "Chaerephon leucogaster"="",
                       "Momopterus acetabulosus"="",
                       "Taphosius mauritianus"="",
                       "Chaerephon sp"="",
                       "Chaerephon pusillus"="",
                       "Pteropus seychellensis comorensis"="",
                       "Rhinolophus lobatus"="",
                       "Rhinolophus rhodesiae"="",
                       "Rhinolophus sp"="",
                       "Pteropus seychellensis seychellensis"="",
                       "Pipistrellus javanicas"="",
                       "Scotphilus kuhlii"="",
                       "Rousettus leschenaulti"="",
                       "Miniopterus schreibersi"="",
                       "Rhinolophus pearsoni"="",
                       "Myotis ricketti"="",
                       "Rhinolophus blythi"="",
                       "Rhinolophus monoceros"="",
                       "Rhinolophus shamelli"="",
                       "Rhinolophus spp"="",
                       "myotis bechsteinii"="",
                       "Neoromicia cf. zuluensis"="",
                       ""="",
                       "Myotis oxygnathus"="",
                       "Myotis mistacinus"="",
                       "Hypsugo savii"="",
                       "Chaerephon pumila"="",
                       "Myotis bocagei"="",
                       "Pipistrellus capensis"="",
                       "Lyssonycteris angolensis"="",
                       "Nycteris argae"="",
                       "Artibeus phaeotis"="",
                       "Artibeus watsoni"="",
                       "Mormoops megalohyla"="",
                       "Pteropus medius"="",
                       "Myonicteris angolensis"="",
                       "Stenonycteris lanosus"="",
                       "Pteropus sp"="",
                       "Microchiroptera (order: chiroptera)"="",
                       "Hipposideros terasensis"="",
                       "Scotorepens rueppellii"="",
                       "Cynopterus spp"="",
                       "Eonycteris spp"="",
                       "Macroglossus spp"="",
                       "Rousettus spp"="",
                       "Taphozous spp"="",
                       "Chalinolobus spp"="",
                       "Scotophilus spp"="",
                       "Scotorepens spp"="",
                       "Chaerephon spp"="",
                       "Molossidae bats"="",
                       "Neuromicia nana"="",
                       "Neuromicia helios"="",
                       "Neoromicia spp"="",
                       "Nycticeinops schlieffenii"="",
                       "Rhinolophus darlingi damarensis"="",
                       "Insectivorous bat species"="",
                       "Rousettous aegyptiacus"="",
                       "Taphozous perforates"="",
                       "Pipistrellus deserti,Rousettous aegyptiacus,Nycteris thebaica,Taphozous perforates"="",
                       "?"="",
                       "Rhinolophus ferrumequinum,Myotis macrodactylus"="",
                       "Myotis aurascens Kuzyakin,Myotis petax"="",
                       "Coelops frithii formosanus"="",
                       "Hipposideros armiger terasensis"="",
                       "Barbastella darjelingesis"="",
                       "Myotis fimbriatus taiwanensis"="",
                       "Myotis rufoniger"="",
                       "Pipistrellus montanus"="",
                       "Pipistrellus taiwanesis"="",
                       "Myotis formosus flavus"="",
                       "Myotis blythii,Hypsugo savii"="",
                       "Fruit bat"="",
                       "(multiple) species of horseshoe bat"="",
                       "Rhinolophus sinicus,Rhinolophus ferrumequinum"="",
                       "Aselliscus stoliczkanus,Rhinolophus affinis"="",
                       "Chaerephon spp."="",
                       "Pipistrellus nanus"="",
                       "Hipposideros caffer ruber"="",
                       "Myotis sp."="",
                       "Rhinolophus rouxi"="",
                       "Pipistrellus sp."="",
                       "Eonycteris spelaea,Rousettus spp,Rousettus leschenaultia"="",
                       "Vespertilio superans"="",
                       "Rhinolophus affiinus"="",
                       "Rousettus lechenaulti"="",
                       "Tyloncyteris pachypus"="",
                       "Myotis emarginatus,Rhinolophus ferrumequinum"="",
                       "Hypsugo pulveratus"="",
                       "Miniopterus filiginosus"="",
                       "Peking Myotis,Myotis pequinius"="",
                       "Rhinolophus malyanus"="",
                       "Coelops frithi"="",
                       "Craseonycteris thonlongyal"="",
                       "Rhinoloplus malayanus"="",
                       "Rhinoloplus sp."="",
                       "Miniopterus meghrebensis"="",
                       "Rhinolophus borneenis"="",
                       "Tadarida spp"="",
                       "Myotis spp"="",
                       "Pipistrellus minus"="",
                       "Pipistrellus spp"="",
                       "Tylonycteris spp"="",
                       "Unclassified"="",
                       "Scotophilus"="",
                       "Cynopterus,Macroglossus"="",
                       "Cynopterus"="",
                       "Cynopterus,Megaerops"="",
                       "Hipposideros,Rhinolophus"="",
                       "Megaderma,Cynopterus,Rhinolophus"="",
                       "Cynopterus,Rhinolophus"="",
                       "Rhinolophus"="",
                       "Pipistrellus"="",
                       "Myotis"="",
                       "Taphozous,Scotophilus"="",
                       "Megaerops,Eonycteris,Rousettus,Rhinolophus"="",
                       "Cynopterus,Pteropus,Rousettus"="",
                       "Cynopterus,Macroglossus,Megaerops"="",
                       "Cynopterus,Tylonycteris"="",
                       "Taphozous,Cynopterus,Megaderma"="",
                       "Taphozous,Hipposideros,Megaderma,Rousettus,Rhinolophus"="",
                       "Cynopterus,Megaerops,Rousettus"="",
                       "Cynopterus,Eonycteris"="",
                       "Cynopterus,Eonycteris,Megaerops"="",
                       "Megaerops"="",
                       "Ia"="",
                       "Hipposideros"="",
                       "Hipposideros,Cynopterus,Eonycteris,Macroglossus,Megaerops,Rousettus"="",
                       "Rousettus"="",
                       "Hipposideros,Cynopterus,Eonycteris,Macroglossus,Megaerops,Rousettus,Rhinolophus"="",
                       "Hipposideros,Megaerops,Rousettus,Rhinolophus"="",
                       "Hipposideros,Eonycteris,Cynopterus,Megaerops,Rousettus,Rhinolophus"="",
                       "Aselliscus,Rhinolophus"="",
                       "Megaerops,Rousettus"="",
                       "Megaerops,Rousettus,Eonycteris"="",
                       "Hipposideros cf. gigas"="",
                       "Hipposideros cf. ruber"="",
                       "Nycteris cf. gambiensis"="",
                       "Lissonycteris angolensis"="",
                       "Murina spp"="",
                       "Megaerops kusnotei"="",
                       "Rhinolophus affinus"="",
                       "Molossus major"="",
                       "Mormoops sp."="",
                       "Molossus molossus,Tadarida brasiliensis"="",
                       "Dobsonia mollucensis"="",
                       "Haplonicteris fischeri"="",
                       "Miniopterus orianae bassanii"="",
                       "Miniopterus orianae oceanensis"="",
                       "Nyctiphilus geoffroyi"="",
                       "Nyctophilus major"="",
                       "Austronomus australis"="",
                       "Ozimops sp"="",
                       "Hipposideros sp"="",
                       "Nycteris sp."="",
                       "Epomps buettikoferi"="",
                       "Myotis oxygnatus"="",
                       "Myotis capaccini"="",
                       "Rhinolophus cornutus"="",
                       "Pipistrelus abramus"="",
                       "Myotis formosus chofukusei"="",
                       "Pteropus spp."="",
                       "Taphozous"="",
                       "Miniopterus africanus"="",
                       "Vampyriscus nymphaea"="",
                       "Glossophaga comissarisi"="",
                       "Eptesicus furnalis"="",
                       "Miniopterus sp."="",
                       "Scotoecus sp."=""))

## check
setdiff(data$species,tree$tip)
length(setdiff(data$species,tree$tip))
length(unique(data$species))

#trim dataset to only phylogeny 
set=data[data$species%in%tree$tip.label,,]
set$species=factor(set$species)
set$studies=factor(set$Field.25)