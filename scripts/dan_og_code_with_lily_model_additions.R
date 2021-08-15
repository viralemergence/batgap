##bat coronavirus gap analysis

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

##set working directory 
setwd("~/Desktop/batgap") 

##load data (in github: batgap/data/batdata.csv)
data <- read.csv("batdata.csv")
which(data$bat_species=="")
data <- data[-c(98,381,382),,]
data$bat_species <- as.character(data$bat_species)

##load in Upham phylogeny from github: batgap/phylos/
tree <- read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")

## load in taxonomy from github: batgap/phylos/
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

## are all bats in tree
setdiff(data$bat_species,tree$tip)

## multiple species?
data$spec=sapply(strsplit(data$bat_species,","),function(x) length(x))

## fix others
data$species=data$bat_species
data$species=revalue(data$species,
                     c("Eptesicus"="drop",
                       "Artibeus literatus"="Artibeus lituratus",
                       "Vampyressa pussilla"="Vampyressa pusilla",
                       "Myotis cillolabrum"="Myotis ciliolabrum",
                       "Macronycteris commersoni"="Hipposideros commersoni",
                       "Chaerephon leucogaster"="Chaerephon pumilus",
                       "Momopterus acetabulosus"="Mormopterus acetabulosus",
                       "Taphosius mauritianus"="Taphozous mauritianus",
                       "Chaerephon sp"="drop",
                       "Chaerephon pusillus"="Chaerephon pumilus",
                       "Pteropus seychellensis comorensis"="Pteropus seychellensis",
                       "Rhinolophus lobatus"="landeri",
                       "Rhinolophus rhodesiae"="Rhinolophus simulator",
                       "Rhinolophus sp"="drop",
                       "Pteropus seychellensis seychellensis"="Pteropus seychellensis",
                       "Pipistrellus javanicas"="Pipistrellus javanicus",
                       "Scotphilus kuhlii"="Scotophilus kuhlii",
                       "Rousettus leschenaulti"="Rousettus leschenaultii",
                       "Miniopterus schreibersi"="Miniopterus schreibersii",
                       "Rhinolophus pearsoni"="Rhinolophus pearsonii",
                       "Myotis ricketti"="Myotis pilosus",
                       "Rhinolophus blythi"="Rhinolophus lepidus",
                       "Rhinolophus monoceros"="Rhinolophus pusillus",
                       "Rhinolophus shamelli"="Rhinolophus shameli",
                       "Rhinolophus spp"="drop",
                       "myotis bechsteinii"="Myotis bechsteinii",
                       "Neoromicia cf. zuluensis"="Neoromicia zuluensis",
                       "Myotis oxygnathus"="Myotis blythii",
                       "Myotis mistacinus"="Myotis mystacinus",
                       "Hypsugo savii"="Pipistrellus savii",
                       "Chaerephon pumila"="Chaerephon pumilus",
                       "Myotis bocagei"="Myotis bocagii",
                       "Pipistrellus capensis"="Neoromicia capensis",
                       "Lyssonycteris angolensis"="Myonycteris angolensis",
                       "Nycteris argae"="Nycteris arge",
                       "Artibeus phaeotis"="Dermanura phaeotis",
                       "Artibeus watsoni"="Dermanura watsoni",
                       "Mormoops megalohyla"="Mormoops megalophylla",
                       "Pteropus medius"="Pteropus giganteus",
                       "Myonicteris angolensis"="Myonycteris angolensis",
                       "Stenonycteris lanosus"="Rousettus lanosus",
                       "Pteropus sp"="drop",
                       "Microchiroptera (order: chiroptera)"="drop",
                       "Hipposideros terasensis"="Hipposideros armiger",
                       "Scotorepens rueppellii"="Scoteanax rueppellii",
                       "Cynopterus spp"="drop",
                       "Eonycteris spp"="drop",
                       "Macroglossus spp"="drop",
                       "Rousettus spp"="drop",
                       "Taphozous spp"="drop",
                       "Chalinolobus spp"="drop",
                       "Scotophilus spp"="drop",
                       "Scotorepens spp"="drop",
                       "Chaerephon spp"="drop",
                       "Molossidae bats"="drop",
                       "Neuromicia nana"="Neoromicia nana",
                       "Neuromicia helios"="Neoromicia helios",
                       "Neoromicia spp"="drop",
                       "Nycticeinops schlieffenii"="Nycticeinops schlieffeni",
                       "Rhinolophus darlingi damarensis"="Rhinolophus darlingi",
                       "Insectivorous bat species"="drop",
                       "Rousettous aegyptiacus"="Rousettus aegyptiacus",
                       "Taphozous perforates"="Taphozous perforatus",
                       "Pipistrellus deserti,Rousettous aegyptiacus,Nycteris thebaica,Taphozous perforates"="drop",
                       "?"="drop",
                       "Rhinolophus ferrumequinum,Myotis macrodactylus"="drop",
                       "Myotis aurascens Kuzyakin,Myotis petax"="drop",
                       "Coelops frithii formosanus"="Coelops frithii",
                       "Hipposideros armiger terasensis"="Hipposideros armiger",
                       "Barbastella darjelingesis"="Barbastella beijingensis",
                       "Myotis fimbriatus taiwanensis"="Myotis fimbriatus",
                       "Myotis rufoniger"="drop",
                       "Pipistrellus montanus"="drop",
                       "Pipistrellus taiwanesis"="drop",
                       "Myotis formosus flavus"="Myotis formosus",
                       "Myotis blythii,Hypsugo savii"="drop",
                       "Fruit bat"="drop",
                       "(multiple) species of horseshoe bat"="drop",
                       "Rhinolophus sinicus,Rhinolophus ferrumequinum"="drop",
                       "Aselliscus stoliczkanus,Rhinolophus affinis"="drop",
                       "Chaerephon spp."="drop",
                       "Pipistrellus nanus"="Neoromicia nana",
                       "Hipposideros caffer ruber"="Hipposideros ruber",
                       "Myotis sp."="drop",
                       "Rhinolophus rouxi"="Rhinolophus rouxii",
                       "Pipistrellus sp."="drop",
                       "Eonycteris spelaea,Rousettus spp,Rousettus leschenaultia"="drop",
                       "Vespertilio superans"="Vespertilio sinensis",
                       "Rhinolophus affiinus"="Rhinolophus affinis",
                       "Rousettus lechenaulti"="Rousettus leschenaultii",
                       "Tyloncyteris pachypus"="Tylonycteris pachypus",
                       "Myotis emarginatus,Rhinolophus ferrumequinum"="drop",
                       "Hypsugo pulveratus"="Pipistrellus pulveratus",
                       "Miniopterus filiginosus"="Miniopterus fuliginosus",
                       "Peking Myotis,Myotis pequinius"="Myotis pequinius",
                       "Rhinolophus malyanus"="Rhinolophus malayanus",
                       "Coelops frithi"="Coelops frithii",
                       "Craseonycteris thonlongyal"="Craseonycteris thonglongyai",
                       "Rhinoloplus malayanus"="Rhinolophus malayanus",
                       "Rhinoloplus sp."="drop",
                       "Miniopterus meghrebensis"="Miniopterus maghrebensis",
                       "Rhinolophus borneenis"="Rhinolophus borneensis",
                       "Tadarida spp"="drop",
                       "Myotis spp"="drop",
                       "Pipistrellus minus"="Pipistrellus tenuis",
                       "Pipistrellus spp"="drop",
                       "Tylonycteris spp"="drop",
                       "Unclassified"="drop",
                       "Scotophilus"="drop",
                       "Cynopterus,Macroglossus"="drop",
                       "Cynopterus"="drop",
                       "Cynopterus,Megaerops"="drop",
                       "Hipposideros,Rhinolophus"="drop",
                       "Megaderma,Cynopterus,Rhinolophus"="drop",
                       "Cynopterus,Rhinolophus"="drop",
                       "Rhinolophus"="drop",
                       "Pipistrellus"="drop",
                       "Myotis"="drop",
                       "Taphozous,Scotophilus"="drop",
                       "Megaerops,Eonycteris,Rousettus,Rhinolophus"="drop",
                       "Cynopterus,Pteropus,Rousettus"="drop",
                       "Cynopterus,Macroglossus,Megaerops"="drop",
                       "Cynopterus,Tylonycteris"="drop",
                       "Taphozous,Cynopterus,Megaderma"="drop",
                       "Taphozous,Hipposideros,Megaderma,Rousettus,Rhinolophus"="drop",
                       "Cynopterus,Megaerops,Rousettus"="drop",
                       "Cynopterus,Eonycteris"="drop",
                       "Cynopterus,Eonycteris,Megaerops"="drop",
                       "Megaerops"="drop",
                       "Ia"="drop",
                       "Hipposideros"="drop",
                       "Hipposideros,Cynopterus,Eonycteris,Macroglossus,Megaerops,Rousettus"="drop",
                       "Rousettus"="drop",
                       "Hipposideros,Cynopterus,Eonycteris,Macroglossus,Megaerops,Rousettus,Rhinolophus"="drop",
                       "Hipposideros,Megaerops,Rousettus,Rhinolophus"="drop",
                       "Hipposideros,Eonycteris,Cynopterus,Megaerops,Rousettus,Rhinolophus"="drop",
                       "Aselliscus,Rhinolophus"="drop",
                       "Megaerops,Rousettus"="drop",
                       "Megaerops,Rousettus,Eonycteris"="drop",
                       "Hipposideros cf. gigas"="Hipposideros gigas",
                       "Hipposideros cf. ruber"="Hipposideros ruber",
                       "Nycteris cf. gambiensis"="Nycteris gambiensis",
                       "Lissonycteris angolensis"="Myonycteris angolensis",
                       "Murina spp"="drop",
                       "Megaerops kusnotei"="Megaerops kusnotoi",
                       "Rhinolophus affinus"="Rhinolophus affinis",
                       "Molossus major"="Molossus molossus",
                       "Mormoops sp."="drop",
                       "Molossus molossus,Tadarida brasiliensis"="drop",
                       "Dobsonia mollucensis"="Dobsonia moluccensis",
                       "Haplonicteris fischeri"="Haplonycteris fischeri",
                       "Miniopterus orianae bassanii"="Miniopterus orianae",
                       "Miniopterus orianae oceanensis"="Miniopterus orianae",
                       "Nyctiphilus geoffroyi"="Nyctophilus geoffroyi",
                       "Nyctophilus major"="Nyctophilus timoriensis",
                       "Austronomus australis"="Tadarida australis",
                       "Ozimops sp"="drop",
                       "Hipposideros sp"="drop",
                       "Nycteris sp."="drop",
                       "Epomps buettikoferi"="Epomops buettikoferi",
                       "Myotis oxygnatus"="Myotis blythii",
                       "Myotis capaccini"="Myotis capaccinii",
                       "Rhinolophus cornutus"="drop",
                       "Pipistrelus abramus"="Pipistrellus abramus",
                       "Myotis formosus chofukusei"="Myotis formosus",
                       "Pteropus spp."="drop",
                       "Taphozous"="drop",
                       "Miniopterus africanus"="Miniopterus inflatus",
                       "Vampyriscus nymphaea"="Vampyressa nymphaea",
                       "Glossophaga comissarisi"="Glossophaga commissarisi",
                       "Eptesicus furnalis"="Eptesicus furinalis",
                       "Miniopterus sp."="drop",
                       "Scotoecus sp."="drop"))

## check
setdiff(data$species,tree$tip)
length(setdiff(data$species,tree$tip))
length(unique(data$species))

#trim dataset to only phylogeny 
set=data[data$species%in%tree$tip.label,,]
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

##add column for bat family 
set$species_char <- as.character(set$species)

which(taxa$species==set$species_char[1])

data_fam <- data.frame(family = taxa$fam[149])

for(i in 2:2272) {
  new_row <- data.frame(family = taxa$fam[which(taxa$species==set$species_char[i])])
  data_fam <- rbind(data_fam, new_row)
}

set <- cbind(set, data_fam)

##add column for country hemisphere (not put in model, but helped me correct sampling seasons)
set$country <- as.character(set$country)
set$country <- revalue(set$country, c("?"="NA"))
set$country <- as.factor(set$country)

set$countryhemi <- set$country 
set$countryhemi <- as.character(set$countryhemi)
set$countryhemi <- revalue(set$countryhemi, c("Brazil"="S","Canada"="N","Italy"="N","United Kingdom"="N",
                                              "France"="N","United States of America"="N","Germany"="N",
                                              "Madagascar"="S", "Mauritius"="S","Mozambique"="S",
                                              "Seychelles"="S","Philippines"="N","China"="N","Netherlands"="N",
                                              "South Africa"="S","Democratic Republic of the Congo"="S",
                                              "Mexico"="N","Sri Lanka"="N","Rwanda"="S","Vietnam"="N",
                                              "Australia"="S","Taiwan"="N","Indonesia"="S","Papua New Guinea"="S",
                                              "East Timor"="S","Malaysia"="N","Egypt"="N","Denmark"="N","Korea"="N",
                                              "Ghana"="N","Singapore"="N","Slovenia"="N","Thailand"="N",
                                              "Saudi Arabia"="N","Kenya"="N","Hungary"="N","New Zealand"="S",
                                              "Belgium"="N","Luxembourg"="N","Spain"="N","Morocco"="N",
                                              "Tunisia"="N","Malaysian Borneo"="N","Finland"="N","Lebanon"="N",
                                              "Trinidad"="N","Gabon"="N","Japan"="N","Nigeria"="N","Guinea"="N",
                                              "Bulgaria"="N","South Korea"="N","Ukraine"="N","Romania"="N",
                                              "Myanmar"="N","Costa Rica"="N"))

#fix sampling season column
#in github: batgap/data/dataforcorrectsampleseason.xlsx
dataforsampleseason <- read_excel("~/Desktop/batgap/dataforcorrectsampleseason.xlsx")
set$sample_season_corrected <- dataforsampleseason$sample_season_corrected
unique(set$sample_season_corrected)
as.character(unique(set$sample_season_corrected))

set$summer <- set$sample_season_corrected
set$summer <- revalue(set$summer, c("Winter"=0,"Summer"=1,"Fall"=0,"Year"=1,"Fall,Winter"=0,"Spring"=0,
                                    "Summer,Fall,Spring"=1,"Summer,Fall"=1,"Winter,Spring"=0,"Spring,Summer,Fall"=1,
                                    "Spring,Summer"=1,"Spring,Summer,Winter"=1,"Winter,Spring,Summer"=1,"Spring,Fall"=0,
                                    "Summer,Spring"=1))
set$fall <- set$sample_season_corrected
set$fall <- revalue(set$fall, c("Winter"=0,"Summer"=0,"Fall"=1,"Year"=1,"Fall,Winter"=1,"Spring"=0,
                                "Summer,Fall,Spring"=1,"Summer,Fall"=1,"Winter,Spring"=0,"Spring,Summer,Fall"=1,
                                "Spring,Summer"=0,"Spring,Summer,Winter"=0,"Winter,Spring,Summer"=0,"Spring,Fall"=1,
                                "Summer,Spring"=0))
set$winter <- set$sample_season_corrected
set$winter <- revalue(set$winter, c("Winter"=1,"Summer"=0,"Fall"=0,"Year"=1,"Fall,Winter"=1,"Spring"=0,
                                    "Summer,Fall,Spring"=0,"Summer,Fall"=0,"Winter,Spring"=1,"Spring,Summer,Fall"=0,
                                    "Spring,Summer"=0,"Spring,Summer,Winter"=1,"Winter,Spring,Summer"=1,"Spring,Fall"=0,
                                    "Summer,Spring"=0))
set$spring <- set$sample_season_corrected
set$spring <- revalue(set$spring, c("Winter"=0,"Summer"=0,"Fall"=0,"Year"=1,"Fall,Winter"=0,"Spring"=1,
                                    "Summer,Fall,Spring"=1,"Summer,Fall"=0,"Winter,Spring"=1,"Spring,Summer,Fall"=1,
                                    "Spring,Summer"=1,"Spring,Summer,Winter"=1,"Winter,Spring,Summer"=1,"Spring,Fall"=1,
                                    "Summer,Spring"=1))

##fix latitudes + longitudes 
angle2dec <- function(angle) {
  angle <- as.character(angle)
  angle <- ifelse(grepl("S|W", angle), paste0("-", angle), angle)
  angle <- trimws(gsub("[^- +.0-9]", "", angle))
  x <- do.call(rbind, strsplit(angle, split=' '))
  x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    (abs(y[1]) + y[2]/60 + y[3]/3600) * sign(y[1])
  })
  return(x)
}

set$latitude <- as.character(set$latitude) 
set$latitude[which(set$latitude=="")] <- "NA"
set$latitude <- revalue(set$latitude, c("51º77'27\"N"=angle2dec("51º 77' 27'' N"),"51º39’96\"N"=angle2dec("51º 39’ 96'' N"),"47º 00'39” N"=angle2dec("47º 00' 39'' N"),"47º 04'59” N"=angle2dec("47º 04' 59'' N"),"46º 56'23” N"=angle2dec("46º 56' 23'' N"),"47º 02'39” N"=angle2dec("47º 02' 39'' N"),"47º 04'24” N"=angle2dec("47º 04' 24'' N"),"46º 51'27” N"=angle2dec("46º 51' 27'' N"),"46º 49'28” N"=angle2dec("46º 49' 28'' N"),
                                        "50°25′46.91′′N"=angle2dec("50° 25' 46.91'' N"),"24°25'35″N"=angle2dec("24° 25' 35'' N"),"24°25'35″N (Miaoli), 23°19'04″N (Dongshan)"="NA","24°25'35″N (Miaoli), 23°19'04″N (Dongshan), 23°48'24″N (Dili)"="NA","24°25'35″N (Miaoli), 23°19'04″N (Dongshan), 23°48'24″N (Dili), 25°11'24″N (Guihui)"="NA","24°31' 55″N (Dong'ao), 24°31'9″N (Nan'ao), 24°25'35″N (Miaoli)"="NA","24°31'55″N (Dong'ao), 24°31'9″N (Nan'ao), 24°25'35″N (Miaoli)"="NA","4°31'9″N (Nan'ao), 24°25'35″N (Miaoli)"="NA","24°51'50″N (Wulai), 23°19'04″N (Dongshan)"="NA","23°50'51″N (Jhutang), 23°34'05″N (Beigang)"="NA",
                                        "23°37'14″N"=angle2dec("23° 37' 14'' N"),"6°42'2.0\"N"=angle2dec("6° 42' 2.0'' N"),"6°41'6.4\"N"=angle2dec("6° 41' 6.4'' N"),"6°32'22.3\"N"=angle2dec("6° 32' 22.3'' N"),"6°58\"N"=angle2dec("6° 0' 58'' N"),"7°43'24.9\"N"=angle2dec("7° 43' 24.9'' N"),"7°43'25.7\"N"=angle2dec("7° 43' 25.7'' N"),"13830018.9\"N"="13.830018","4391037 N"="39.639155","46°47′S"=angle2dec("46° 47′ 0'' S"),
                                        "54°06'04\"N"=angle2dec("54° 06' 04'' N"),"54°16'44\"N"=angle2dec("54° 16' 44'' N"),"53°56′09′′N"=angle2dec("53° 56' 09'' N"),"54°16'18\"N"=angle2dec("54° 16' 18'' N"),"1°07N"=angle2dec("1° 07' 0'' N"),"0°82N"=angle2dec("0° 82' 0'' N"),"1°36 S"=angle2dec("1° 36' 0'' S"),"0°98N"=angle2dec("0° 98' 0'' N"),"0°86 S"=angle2dec("0° 86' 0'' S"),"23º50'51\"N"=angle2dec("23º 50' 51'' N"),
                                        "23º34'05\"N"=angle2dec("23º 34' 05'' N"),"2016"="NA","42°5'30.7\"N"=angle2dec("42° 5' 30.7'' N"),"42°4'N"=angle2dec("42° 4' 0'' N"),"42°0'21.0\"N"=angle2dec("42° 0' 21.0'' N"),"42°9'32.0\"N"=angle2dec("42° 9' 32.0'' N"),"42°9'7.0\"N"=angle2dec("42° 9' 7.0'' N"),"6°41′13.56′′ N"=angle2dec("6° 41' 13.56'' N"),"7°43′24.899′′ N"=angle2dec("7° 43' 24.899'' N"),"7°35′23.1′′ N"=angle2dec("7° 35' 23.1'' N"),
                                        "6°58′0.001′′ N"=angle2dec("6° 58' 0.001'' N"),"7°15′43.099′′ N"=angle2dec("7° 15' 43.099'' N"),"5°55′44.4′′ N"=angle2dec("5° 55' 44.4'' N"),"50°27′0.324′′ N"=angle2dec("50° 27' 0.324'' N"),"54°14′51.271′′ N"=angle2dec("54° 14' 51.271'' N"),"45°12′0.00′′ N"=angle2dec("45° 12' 0.00'' N"),"50°20′5.316′′ N"=angle2dec("50° 20' 5.316'' N"),"52°1′46.859′′ N"=angle2dec("52° 1' 46.859'' N"),"10º21'14.2\"N, 10º50'08.5\"N"="NA","10º21'14.2\"N"=angle2dec("10º 21' 14.2'' N"),
                                        "10º21'14.2\"N, 10º03'35.6\"W"="NA","9º3226.0\"N"=angle2dec("9º 0' 3226.0'' N"),"10º25'01.4\" N, 9º3226.0\"N"="NA","10º17'31.9\"N"=angle2dec("10º 17' 31.9'' N"),"9º3226.0\"N, 10º50'08.5\"N"="NA","9º56'19.3\"N"=angle2dec("9º 56' 19.3'' N"),"10º17'31.9\"N, 9º56'19.3\"N"="NA","10º21'14.2\"N, 10º17'31.9\"N, 9º56'19.3\"N"="NA","10º25'01.4\" N"=angle2dec("10º 25' 01.4'' N"),"10º21'05.6\" N"=angle2dec("10º 21' 05.6'' N"),
                                        "10º25'01.4\" N, 10º24'11.0\"N"="NA","10º24'11.0\"N"=angle2dec("10º 24' 11.0'' N"),"10º24'11.0\"N, 10º8'35.64\"N"="NA","9º58'39.4\" N"=angle2dec("9º 58' 39.4'' N")))

set$latitude[which(set$latitude=="24°31'9″N (Nan'ao), 24°25'35″N (Miaoli)")] <- "NA"
set$latitude[which(set$latitude=="13830018.9\"N")] <- "13.830018"
set$latitude <- as.numeric(set$latitude)

set$longitude <- as.character(set$longitude)
set$longitude[which(set$longitude=="")] <- "NA"
set$longitude <- revalue(set$longitude, c("1º33'41\"E"=angle2dec("1º 33' 41'' E"),"1º67'75\"E"=angle2dec("1º 67' 75'' E"),"02º 34'54'' E"=angle2dec("02º 34' 54'' E")," 02º 31'02” E"=angle2dec("02º 31' 02'' E"),"02º 32'19” E"=angle2dec("02º 32' 19'' E"),"02º 33'38” E"=angle2dec("02º 33' 38'' E"),"02º 31122” E"=angle2dec("02º 31' 22'' E"),"02º 18'59” E"=angle2dec("02º 18' 59” E"),"02º 23'21” E"=angle2dec("02º 23' 21'' E"),
                                          "6°55′52.17′′E"=angle2dec("6° 55' 52.17'' E"),"121°00'45″E"=angle2dec("121° 00' 45'' E"),"121°00'45″E (Miaoli), 120°25'27″E (Dongshan)"="NA","121°00'45″E (Miaoli), 120°250 27″E (Dongshan), 120°54'51″E (Dili)"="NA","121°00'45″E (Miaoli), 120°25'27″E (Dongshan), 120°54'51″E (Dili), 121°40'41″E (Guihui)"="NA","121°49'53″E (Dong'ao), 121°46'2″E (Nan'ao), 121°00'45″E (Miaoli)"="NA","121°46'2″E (Nan'ao), 121°00'45″E (Miaoli)"="NA","121°33'05″E (Wulai), 120°25'27″E (Dongshan)"="NA","120°23'12″E"=angle2dec("120° 23' 12'' E"),"120°23'12″E (Jhutang), 120°17'51″E (Beigang)"="NA",
                                          "120°15'40″E"=angle2dec("120° 15' 40'' E"),"1°37\"29.9\"W"=angle2dec("1° 37' 29.9'' W"),"1°33'42.8\"W"=angle2dec("1° 33' 42.8'' W"),"1°24'41.5\"W"=angle2dec("1° 24' 41.5'' W"),"1°16\"W"=angle2dec("1° 0' 16'' W")," 1°59'16.5\"W"=angle2dec("1° 59' 16.5'' W"),"1°59'33.5\"W"=angle2dec("1° 59' 33.5'' W"),"101809054.9\"E"="101.8090549","726475 E"="-72.360832","167°38′E"=angle2dec("167° 38' 0'' E"),
                                          "10°14'14\"E"=angle2dec("10° 14' 14'' E"),"9°58'17\"E"=angle2dec("9° 58' 17'' E"),"10°18′57′′E"=angle2dec("10° 18' 57'' E"),"10°17'14\"E"=angle2dec("10° 17' 14'' E"),"13°20 E"=angle2dec("13° 20' 0'' E"),"13°45 E"=angle2dec("13° 45' 0'' E"),"13°46 E"=angle2dec("13° 46' 0'' E")," 13°19 E"=angle2dec("13° 19' 0'' E"),"12°77 E"=angle2dec("12° 77' 0'' E"),"120º23'12\"E"=angle2dec("120º 23' 12'' E"),
                                          "120º17'51\"E"=angle2dec("120º 17' 51'' E"),"27°12'31.3\"E\005"=angle2dec("27° 12' 31.3'' E"),"27°46'E"=angle2dec("27° 46' 0'' E")," 27°25'21.3\"E"=angle2dec("27° 25' 21.3'' E"),"27°30'45.0\"E"=angle2dec("27° 30' 45.0'' E"),"27°21'28.0\"E"=angle2dec("27° 21' 28.0'' E"),"1°20′38.94′′ W"=angle2dec("1° 20' 38.94'' W"),"1°59′16.501′′ W"=angle2dec("1° 59' 16.501'' W"),"1°52′30.299′′ W"=angle2dec("1° 52' 30.299'' W"),"1°16′0.001\" W"=angle2dec("1° 16' 0.001'' W"),
                                          "0°29′29.501\" E"=angle2dec("0° 29' 29.501'' E"),"0°4′30′′ E"=angle2dec("0° 4' 30'' E"),"30°31′24.24\" E"=angle2dec("30° 31' 24.24'' E"),"10°4′3.347′′ E"=angle2dec("10° 4' 3.347'' E"),"29°0′0.00′′ E"=angle2dec("29° 0' 0.00'' E"),"7°14′30.912′′ E"=angle2dec("7° 14' 30.912'' E"),"6°13′4.908′′ E"=angle2dec("6° 13' 4.908'' E"),"85º19'45.7\"W, 85º37'27.2\" W"="NA","85º19'45.7\"W"=angle2dec("85º 19' 45.7'' W"),"85º19'45.7\"W, 85º15'55.2\"W"="NA",
                                          "84º24'02.5\"W"=angle2dec("84º 24' 02.5'' W"),"84º07'25.4\"W, 84º24'02.5\"W"="NA","84º07'21.3\"W"=angle2dec("84º 07' 21.3'' W"),"84º24'02.5\"W, 85º37'27.2\" W"="NA","84º02'42.8\"W"=angle2dec("84º 02' 42.8'' W"),"84º07'21.3\"W, 84º02'42.8\"W"="NA","85º19'45.7\"W, 84º07'21.3\"Wm 84º02'42.8\"W"="NA","84º07'25.4\"W"=angle2dec("84º 07' 25.4'' W"),"85º25'05.8 W"=angle2dec("85º 25' 5.8'' W"),"84º07'25.4\"W, 84º08'01.6\"W"="NA",
                                          "84º08'01.6\"W"=angle2dec("84º 08' 01.6'' W"),"84º08'01.6\"W, 85º26'47.33\" W"="NA","84º50'16.7\"W"=angle2dec("84º 50' 16.7'' W")))
set$longitude <- as.numeric(set$longitude)

###models###

##model with no covariates
model=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix),
             method="REML",mods=~1,data=set,
             control=list(optimizer="optim", optmethod="BFGS"))

summary(model)

##model with study_type, PCR method, tissue_simplified + detection_method as mods 
model_with_mods=rma.mv(yi=yi,V=vi,
                       random=list(~1|study/observation,~1|species,~1|phylo),
                       R=list(phylo=cmatrix),
                       method="REML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method, data=set,
                       control=list(optimizer="optim", optmethod="BFGS"))

summary(model_with_mods)

##use anova to characterize difference between model and model_with_mods
anova(model_with_mods, btt=2)
anova(update(model_with_mods,method="ML"),update(model,method="ML"))

##get R2ish thing (Dan to get citation for function below)
r2fun=function(model,model_with_mods){
  
  ## both model (no moderators) and model with moderators must be fit to (1) the same dataset and (2) with REML
  pr2=(sum(model$sigma2) - sum(model_with_mods$sigma2)) / sum(model$sigma2)  
  pr2=ifelse(pr2<0,0,pr2)
  return(pr2)
}

r2fun(model,model_with_mods)

##more models! 

##model with detection type as mod
model_with_detectiontype=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix),
                                method="REML",mods=~detection_type,data=set,
                                control=list(optimizer="optim", optmethod="BFGS"))
summary(model_with_detectiontype)

##model with bat family + sampling season as mods 
model_fam_season=rma.mv(yi=yi,V=vi,
                        random=list(~1|study/observation,~1|species,~1|phylo),
                        R=list(phylo=cmatrix),
                        method="REML",mods=~family + summer + winter + fall + spring, data=set,
                        control=list(optimizer="optim", optmethod="BFGS"))

summary(model_fam_season)

##model with bat family + sampling season + lat/longs as mods 
model_fam_season_latlong=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix),
                                method="REML",mods=~family + summer + winter + fall + spring + latitude + longitude, data=set,
                                control=list(optimizer="optim", optmethod="BFGS"))
summary(model_fam_season_latlong)

###rerun models above with just alpha + just beta coronaviruses

#####first, new datasets with just alpha + just beta prevalences

set$virus_genus <- as.character(set$virus_genus)

## names of viral genera that aren't strictly alpha or beta
vnames=c("alphacoronavirus or betacoronavirus or gammacoronavirus",
         "alphacoronavirus or betacoronavirus",
         "alphacoronavirus or betacoronavirus or gammacoronavirus or deltacoronavirus",
         "?",
         "alphacoronavirus/betacoronavirus coinfection",
         "alphacoronavirus and betacoronavirus",
         "alphacoronavirus and betacoronavirus and independent bat coronavirus",
         "independent bat coronavirus")

## flag as alphaCoV
set$alpha=ifelse(set$virus_genus%in%c("alphacoronavirus",
                                      "alphacoronavirus/alphacoronavirus coinfection"),1,0)

## count any zeros using general detection as alpha-negative
set$alpha=ifelse(set$positives==0 & set$virus_genus%in%vnames,1,set$alpha)

## repeat with beta
set$beta=ifelse(set$virus_genus%in%c("betacoronavirus","betacoronavirus "),1,0)
set$beta=ifelse(set$positives==0 & set$virus_genus%in%vnames,1,set$beta)

## check
cset=set[c("virus_genus","positives","sample","alpha","beta")]

## make positive columns for alpha and beta
set$positives_alpha=ifelse(set$alpha==1,set$positives,NA)
set$positives_beta=ifelse(set$beta==1,set$positives,NA)

## make sample columns for alpha and beta
set$sample_alpha=ifelse(set$alpha==1,set$sample,NA)
set$sample_beta=ifelse(set$beta==1,set$sample,NA)

## pft in escalc for yi and vi 
set=data.frame(set,escalc(xi=set$positives,ni=set$sample,measure="PFT"))

## back transform
set$backtrans=transf.ipft(set$yi,set$sample)

## check
plot(set$prevalence,set$backtrans)
abline(0,1)

## make alpha and beta yi and vi
alpha_es=escalc(xi=set$positives_alpha,ni=set$sample_alpha,measure="PFT")
beta_es=escalc(xi=set$positives_beta,ni=set$sample_beta,measure="PFT")

## append names
names(alpha_es)=paste(names(alpha_es),"_alpha",sep="")
names(beta_es)=paste(names(beta_es),"_beta",sep="")

## cbind
set=data.frame(set,alpha_es,beta_es)

set_alpha=set[!is.na(set$yi_alpha),]
set_beta=set[!is.na(set$yi_beta),]

######now, rerun models with set_alpha + set_beta 

##model with no covariates: alpha
model_alpha=rma.mv(yi=yi,V=vi,
             random=list(~1|study/observation,~1|species,~1|phylo),
             R=list(phylo=cmatrix),
             method="REML",mods=~1,data=set_alpha,
             control=list(optimizer="optim", optmethod="BFGS"))

summary(model_alpha)

##model with no covariates: beta
model_beta=rma.mv(yi=yi,V=vi,
                   random=list(~1|study/observation,~1|species,~1|phylo),
                   R=list(phylo=cmatrix),
                   method="REML",mods=~1,data=set_beta,
                   control=list(optimizer="optim", optmethod="BFGS"))

summary(model_beta)

##model with study_type, PCR method, tissue_simplified + detection_method as mods: alpha
model_with_mods_alpha=rma.mv(yi=yi,V=vi,
                       random=list(~1|study/observation,~1|species,~1|phylo),
                       R=list(phylo=cmatrix),
                       method="REML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method, data=set_alpha,
                       control=list(optimizer="optim", optmethod="BFGS"))

summary(model_with_mods_alpha)

##model with study_type, PCR method, tissue_simplified + detection_method as mods: beta
model_with_mods_beta=rma.mv(yi=yi,V=vi,
                       random=list(~1|study/observation,~1|species,~1|phylo),
                       R=list(phylo=cmatrix),
                       method="REML",mods=~study_type + single.multiple.PCR.method + tissue_simplified + detection_method, data=set_beta,
                       control=list(optimizer="optim", optmethod="BFGS"))

summary(model_with_mods_beta)

##model with detection type as mod: alpha
model_with_detectiontype_alpha=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix),
                                method="REML",mods=~detection_type,data=set_alpha,
                                control=list(optimizer="optim", optmethod="BFGS"))
summary(model_with_detectiontype_alpha)

##model with detection type as mod: beta
model_with_detectiontype_beta=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix),
                                method="REML",mods=~detection_type,data=set_beta,
                                control=list(optimizer="optim", optmethod="BFGS"))
summary(model_with_detectiontype_beta)

##model with bat family + sampling season as mods: alpha
model_fam_season_alpha=rma.mv(yi=yi,V=vi,
                        random=list(~1|study/observation,~1|species,~1|phylo),
                        R=list(phylo=cmatrix),
                        method="REML",mods=~family + summer + winter + fall + spring, data=set_alpha,
                        control=list(optimizer="optim", optmethod="BFGS"))

summary(model_fam_season_alpha)

##model with bat family + sampling season as mods: beta 
model_fam_season_beta=rma.mv(yi=yi,V=vi,
                        random=list(~1|study/observation,~1|species,~1|phylo),
                        R=list(phylo=cmatrix),
                        method="REML",mods=~family + summer + winter + fall + spring, data=set_beta,
                        control=list(optimizer="optim", optmethod="BFGS"))

summary(model_fam_season_beta)

##model with bat family + sampling season + lat/longs as mods: alpha
model_fam_season_latlong_alpha=rma.mv(yi=yi,V=vi,
                                random=list(~1|study/observation,~1|species,~1|phylo),
                                R=list(phylo=cmatrix),
                                method="REML",mods=~family + summer + winter + fall + spring + latitude + longitude, data=set_alpha,
                                control=list(optimizer="optim", optmethod="BFGS"))

summary(model_fam_season_latlong_alpha)

##model with bat family + sampling season + lat/longs as mods: beta
model_fam_season_latlong_beta=rma.mv(yi=yi,V=vi,
                                      random=list(~1|study/observation,~1|species,~1|phylo),
                                      R=list(phylo=cmatrix),
                                      method="REML",mods=~family + summer + winter + fall + spring + latitude + longitude, data=set_beta,
                                      control=list(optimizer="optim", optmethod="BFGS"))

summary(model_fam_season_latlong_beta)