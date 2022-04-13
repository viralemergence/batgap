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
library(cowplot)

##set working directory 
setwd("~/Documents/GitHub/batgap/data")

##Airtable as of April 13, 2022
data <- read.csv("airtable.csv")

#add column with simplified tissues
unique(data$tissue)
data$tissue_simplified<-data$tissue
data$tissue_simplified<-revalue(data$tissue_simplified,
                                c("rectal swab"="faecal, rectal, or anal",
                                  "oral swab"="oropharyngeal",
                                  "blood"="blood or serum",
                                  "enteric content"="intestine",
                                  "intestines"="intestine",
                                  "lung"="lung or respiratory",
                                  "faeces"="faecal, rectal, or anal",
                                  "liver"="liver",
                                  "faeces/anal swab"="faecal, rectal, or anal",
                                  "rectal/oral swab"="pooled swabs/samples",
                                  "serum"="blood or serum",
                                  "faecal swab"="faecal, rectal, or anal",
                                  "respiratory swab"="lung or respiratory",
                                  "faecal swab/faeces"="faecal, rectal, or anal",
                                  "alimentary specimen"="intestine",
                                  "respiratory specimen"="lung or respiratory",
                                  "faeces or urine swab"="pooled swabs/samples",
                                  "urine swab"="urinary",
                                  "faecal and/or throat sample"="pooled swabs/samples",
                                  "rectal/oral swab,serum"="pooled swabs/samples",
                                  "urine"="urinary",
                                  "faeces or anal swab or intestinal content or oropharyngeal swab"="pooled swabs/samples",
                                  "faeces or rectal tissue/swab"="faecal, rectal, or anal",
                                  "throat swab"="oropharyngeal",
                                  "anal swab"="faecal, rectal, or anal",
                                  "roost feces"="faecal, rectal, or anal",
                                  "roost feces,faeces"="faecal, rectal, or anal",
                                  "pooled tissue (liver/lung/small intestine/brain/kidney/spleen)"="pooled tissue",
                                  "faeces or urine"="pooled swabs/samples",
                                  "fecal swab"="faecal, rectal, or anal",
                                  "fecal pellets"="faecal, rectal, or anal",
                                  "\"saliva, faeces, and urine samples\""="pooled swabs/samples",
                                  "anal swabs/faecal samples"="faecal, rectal, or anal",
                                  "faecal samples"="faecal, rectal, or anal",
                                  "anal swab,nasopharyngeal swabs"="pooled swabs/samples",
                                  "skin swab"="skin swab",
                                  "guano"="faecal, rectal, or anal",
                                  "oral swab,alimentary specimen"="pooled swabs/samples",
                                  "intestines,spleen"="pooled tissue",
                                  "fecal pellets,anal swab"="faecal, rectal, or anal",
                                  "anal swab,oral swab"="pooled swabs/samples",
                                  "faecal swab/pellet"="faecal, rectal, or anal",
                                  "faeces or rectal swab"="faecal, rectal, or anal",
                                  "oral swab,rectal swab,serum"="pooled swabs/samples",
                                  "intestines,lung"="pooled tissue",
                                  "nasopharyngeal swabs,faecal swab"="pooled swabs/samples",
                                  "pharyngeal/anal swab"="pooled swabs/samples",
                                  "Carcass"="pooled tissue",
                                  "oral swab,\"tissue (lung,liver,spleen)\""="pooled tissue",
                                  "nasopharyngeal swabs,anal swab"="pooled swabs/samples",
                                  "fecal pellets,rectal swab,oral swab"="pooled swabs/samples",
                                  "fecal pellets,oral swab"="pooled swabs/samples",
                                  "rectal swab,oral swab"="pooled swabs/samples",
                                  "faeces,oral swab"="pooled swabs/samples",
                                  "oral swab,rectal swab"="pooled swabs/samples",
                                  "respiratory/fecal swab"="pooled swabs/samples",
                                  "oral swab,fecal swab"="pooled swabs/samples"))
unique(data$tissue_simplified)
unique(data$country)
data$country <- revalue(data$country, c("?"="NA"))

#add seasons
data$summer <- data$sample_season_simplified
data$summer <- revalue(data$summer, c("Winter"=0,"Summer"=1,"Fall"=0,"Year"=1,"Fall,Winter"=0,"Spring"=0,
                                    "Summer,Fall,Spring"=1,"Summer,Fall"=1,"Winter,Spring"=0,"Spring,Summer,Fall"=1,
                                    "Spring,Summer"=1,"Spring,Summer,Winter"=1,"Winter,Spring,Summer"=1,"Spring,Fall"=0,
                                    "Summer,Spring"=1))
data$fall <- data$sample_season_simplified
data$fall <- revalue(data$fall, c("Winter"=0,"Summer"=0,"Fall"=1,"Year"=1,"Fall,Winter"=1,"Spring"=0,
                                "Summer,Fall,Spring"=1,"Summer,Fall"=1,"Winter,Spring"=0,"Spring,Summer,Fall"=1,
                                "Spring,Summer"=0,"Spring,Summer,Winter"=0,"Winter,Spring,Summer"=0,"Spring,Fall"=1,
                                "Summer,Spring"=0))
data$winter <- data$sample_season_simplified
data$winter <- revalue(data$winter, c("Winter"=1,"Summer"=0,"Fall"=0,"Year"=1,"Fall,Winter"=1,"Spring"=0,
                                    "Summer,Fall,Spring"=0,"Summer,Fall"=0,"Winter,Spring"=1,"Spring,Summer,Fall"=0,
                                    "Spring,Summer"=0,"Spring,Summer,Winter"=1,"Winter,Spring,Summer"=1,"Spring,Fall"=0,
                                    "Summer,Spring"=0))
data$spring <- data$sample_season_simplified
data$spring <- revalue(data$spring, c("Winter"=0,"Summer"=0,"Fall"=0,"Year"=1,"Fall,Winter"=0,"Spring"=1,
                                    "Summer,Fall,Spring"=1,"Summer,Fall"=0,"Winter,Spring"=1,"Spring,Summer,Fall"=1,
                                    "Spring,Summer"=1,"Spring,Summer,Winter"=1,"Winter,Spring,Summer"=1,"Spring,Fall"=1,
                                    "Summer,Spring"=1))

#fix lat/long
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

data$latitude[which(data$latitude=="")] <- "NA"
data$longitude[which(data$longitude=="")] <- "NA"

unique(data$latitude)
data$latidudedecimal<-as.character(data$latitude)
data$latidudedecimal <- revalue(data$latidudedecimal, c("51º77'27\"N"=angle2dec("51º 77' 27'' N"),"51º39’96\"N"=angle2dec("51º 39’ 96'' N"),"47º 00'39” N"=angle2dec("47º 00' 39'' N"),"47º 04'59” N"=angle2dec("47º 04' 59'' N"),"46º 56'23” N"=angle2dec("46º 56' 23'' N"),"47º 02'39” N"=angle2dec("47º 02' 39'' N"),"47º 04'24” N"=angle2dec("47º 04' 24'' N"),"46º 51'27” N"=angle2dec("46º 51' 27'' N"),"46º 49'28” N"=angle2dec("46º 49' 28'' N"),
                                        "50°25′46.91′′N"=angle2dec("50° 25' 46.91'' N"),"24°25'35″N"=angle2dec("24° 25' 35'' N"),"24°25'35″N (Miaoli), 23°19'04″N (Dongshan)"="NA","24°25'35″N (Miaoli), 23°19'04″N (Dongshan), 23°48'24″N (Dili)"="NA","24°25'35″N (Miaoli), 23°19'04″N (Dongshan), 23°48'24″N (Dili), 25°11'24″N (Guihui)"="NA","24°31' 55″N (Dong'ao), 24°31'9″N (Nan'ao), 24°25'35″N (Miaoli)"="NA","24°31'55″N (Dong'ao), 24°31'9″N (Nan'ao), 24°25'35″N (Miaoli)"="NA","24°31'9″N (Nan'ao), 24°25'35″N (Miaoli)"="NA","24°51'50″N (Wulai), 23°19'04″N (Dongshan)"="NA","24°31' 55″N"=angle2dec("24º 31' 55'' N"),"23°50'51″N (Jhutang), 23°34'05″N (Beigang)"="NA",
                                        "23°37'14″N"=angle2dec("23° 37' 14'' N"),"6°42'2.0\"N"=angle2dec("6° 42' 2.0'' N"),"6°41'6.4\"N"=angle2dec("6° 41' 6.4'' N"),"6°32'22.3\"N"=angle2dec("6° 32' 22.3'' N"),"6°58\"N"=angle2dec("6° 0' 58'' N"),"7°43'24.9\"N"=angle2dec("7° 43' 24.9'' N"),"7°43'25.7\"N"=angle2dec("7° 43' 25.7'' N"),"13830018.9\"N"="13.830018","46°47′S"=angle2dec("46° 47′ 0'' S"),
                                        "54°06'04\"N"=angle2dec("54° 06' 04'' N"),"54°16'44\"N"=angle2dec("54° 16' 44'' N"),"53°56′09′′N"=angle2dec("53° 56' 09'' N"),"54°16'18\"N"=angle2dec("54° 16' 18'' N"),"1°07N"=angle2dec("1° 07' 0'' N"),"0°82N"=angle2dec("0° 82' 0'' N"),"1°36 S"=angle2dec("1° 36' 0'' S"),"0°98N"=angle2dec("0° 98' 0'' N"),"0°86 S"=angle2dec("0° 86' 0'' S"),"23º50'51\"N"=angle2dec("23º 50' 51'' N"),
                                        "23º34'05\"N"=angle2dec("23º 34' 05'' N"),"38.3861ºS"=angle2dec("38.3861° 0' 0'' S"),"37.6515 ºS"=angle2dec("37.6515° 0' 0'' S"),"2016"="NA","42°5'30.7\"N"=angle2dec("42° 5' 30.7'' N"),"42°4'N"=angle2dec("42° 4' 0'' N"),"42°0'21.0\"N"=angle2dec("42° 0' 21.0'' N"),"42°9'32.0\"N"=angle2dec("42° 9' 32.0'' N"),"42°9'7.0\"N"=angle2dec("42° 9' 7.0'' N"),"39.6019\"N (Ogayu), 39.6771\" N (Isari), "="NA","6°41′13.56′′ N"=angle2dec("6° 41' 13.56'' N"),"7°43′24.899′′ N"=angle2dec("7° 43' 24.899'' N"),"7°35′23.1′′ N"=angle2dec("7° 35' 23.1'' N"),
                                        "6°58′0.001′′ N"=angle2dec("6° 58' 0.001'' N"),"7°15′43.099′′ N"=angle2dec("7° 15' 43.099'' N"),"5°55′44.4′′ N"=angle2dec("5° 55' 44.4'' N"),"50°27′0.324′′ N"=angle2dec("50° 27' 0.324'' N"),"54°14′51.271′′ N"=angle2dec("54° 14' 51.271'' N"),"45°12′0.00′′ N"=angle2dec("45° 12' 0.00'' N"),"50°20′5.316′′ N"=angle2dec("50° 20' 5.316'' N"),"52°1′46.859′′ N"=angle2dec("52° 1' 46.859'' N"),"10º21'14.2\"N, 10º50'08.5\"N"="NA","10º21'14.2\"N"=angle2dec("10º 21' 14.2'' N"),
                                        "10º21'14.2\"N, 10º03'35.6\"W"="NA","9º3226.0\"N"=angle2dec("9º 0' 3226.0'' N"),"10º25'01.4\" N, 9º3226.0\"N"="NA","10º17'31.9\"N"=angle2dec("10º 17' 31.9'' N"),"9º3226.0\"N, 10º50'08.5\"N"="NA","9º56'19.3\"N"=angle2dec("9º 56' 19.3'' N"),"10º17'31.9\"N, 9º56'19.3\"N"="NA","10º21'14.2\"N, 10º17'31.9\"N, 9º56'19.3\"N"="NA","10º25'01.4\" N"=angle2dec("10º 25' 01.4'' N"),"10º21'05.6\" N"=angle2dec("10º 21' 05.6'' N"),
                                        "10º25'01.4\" N, 10º24'11.0\"N"="NA","10º24'11.0\"N"=angle2dec("10º 24' 11.0'' N"),"10º24'11.0\"N, 10º8'35.64\"N"="NA","9º58'39.4\" N"=angle2dec("9º 58' 39.4'' N")))
unique(data$latidudedecimal)
data$latidudedecimal <- as.numeric(data$latidudedecimal)

unique(data$longitude)
data$longitudedecimal <- data$longitude
data$longitudedecimal <- revalue(data$longitudedecimal, c("1º33'41\"E"=angle2dec("1º 33' 41'' E"),"1º67'75\"E"=angle2dec("1º 67' 75'' E"),"02º 34'54'' E"=angle2dec("02º 34' 54'' E")," 02º 31'02” E"=angle2dec("02º 31' 02'' E"),"02º 32'19” E"=angle2dec("02º 32' 19'' E"),"02º 33'38” E"=angle2dec("02º 33' 38'' E"),"02º 31122” E"=angle2dec("02º 31' 22'' E"),"02º 18'59” E"=angle2dec("02º 18' 59” E"),"02º 23'21” E"=angle2dec("02º 23' 21'' E"),
                                          "6°55′52.17′′E"=angle2dec("6° 55' 52.17'' E"),"121°49'53″E"=angle2dec("121° 49' 53'' E"),"121°00'45″E"=angle2dec("121° 00' 45'' E"),"121°00'45″E (Miaoli), 120°25'27″E (Dongshan)"="NA","121°00'45″E (Miaoli), 120°250 27″E (Dongshan), 120°54'51″E (Dili)"="NA","121°00'45″E (Miaoli), 120°25'27″E (Dongshan), 120°54'51″E (Dili), 121°40'41″E (Guihui)"="NA","121°49'53″E (Dong'ao), 121°46'2″E (Nan'ao), 121°00'45″E (Miaoli)"="NA","121°46'2″E (Nan'ao), 121°00'45″E (Miaoli)"="NA","121°33'05″E (Wulai), 120°25'27″E (Dongshan)"="NA","120°23'12″E"=angle2dec("120° 23' 12'' E"),"120°23'12″E (Jhutang), 120°17'51″E (Beigang)"="NA",
                                          "120°15'40″E"=angle2dec("120° 15' 40'' E"),"1°37\"29.9\"W"=angle2dec("1° 37' 29.9'' W"),"1°33'42.8\"W"=angle2dec("1° 33' 42.8'' W"),"1°24'41.5\"W"=angle2dec("1° 24' 41.5'' W"),"1°16\"W"=angle2dec("1° 0' 16'' W")," 1°59'16.5\"W"=angle2dec("1° 59' 16.5'' W"),"1°59'33.5\"W"=angle2dec("1° 59' 33.5'' W"),"101809054.9\"E"="101.8090549","167°38′E"=angle2dec("167° 38' 0'' E"),
                                          "10°14'14\"E"=angle2dec("10° 14' 14'' E"),"9°58'17\"E"=angle2dec("9° 58' 17'' E"),"10°18′57′′E"=angle2dec("10° 18' 57'' E"),"10°17'14\"E"=angle2dec("10° 17' 14'' E"),"13°20 E"=angle2dec("13° 20' 0'' E"),"13°45 E"=angle2dec("13° 45' 0'' E"),"13°46 E"=angle2dec("13° 46' 0'' E")," 13°19 E"=angle2dec("13° 19' 0'' E"),"12°77 E"=angle2dec("12° 77' 0'' E"),"120º23'12\"E"=angle2dec("120º 23' 12'' E"),
                                          "120º17'51\"E"=angle2dec("120º 17' 51'' E"),"142.5931ºE"=angle2dec("142.5931° 0' 0'' E"),"145.3173ºE"=angle2dec("145.3173° 0' 0'' E"),"27°12'31.3\"E\005"=angle2dec("27° 12' 31.3'' E"),"27°46'E"=angle2dec("27° 46' 0'' E")," 27°25'21.3\"E"=angle2dec("27° 25' 21.3'' E"),"27°30'45.0\"E"=angle2dec("27° 30' 45.0'' E"),"27°21'28.0\"E"=angle2dec("27° 21' 28.0'' E"),"141.2508\"E (Ogayu),  141.0881\"E (Isari)"='NA',"1°20′38.94′′ W"=angle2dec("1° 20' 38.94'' W"),"1°59′16.501′′ W"=angle2dec("1° 59' 16.501'' W"),"1°52′30.299′′ W"=angle2dec("1° 52' 30.299'' W"),"1°16′0.001\" W"=angle2dec("1° 16' 0.001'' W"),
                                          "0°29′29.501\" E"=angle2dec("0° 29' 29.501'' E"),"0°4′30′′ E"=angle2dec("0° 4' 30'' E"),"30°31′24.24\" E"=angle2dec("30° 31' 24.24'' E"),"10°4′3.347′′ E"=angle2dec("10° 4' 3.347'' E"),"29°0′0.00′′ E"=angle2dec("29° 0' 0.00'' E"),"7°14′30.912′′ E"=angle2dec("7° 14' 30.912'' E"),"6°13′4.908′′ E"=angle2dec("6° 13' 4.908'' E"),"85º19'45.7\"W, 85º37'27.2\" W"="NA","85º19'45.7\"W"=angle2dec("85º 19' 45.7'' W"),"85º19'45.7\"W, 85º15'55.2\"W"="NA",
                                          "84º24'02.5\"W"=angle2dec("84º 24' 02.5'' W"),"84º07'25.4\"W, 84º24'02.5\"W"="NA","84º07'21.3\"W"=angle2dec("84º 07' 21.3'' W"),"84º24'02.5\"W, 85º37'27.2\" W"="NA","84º02'42.8\"W"=angle2dec("84º 02' 42.8'' W"),"84º07'21.3\"W, 84º02'42.8\"W"="NA","85º19'45.7\"W, 84º07'21.3\"Wm 84º02'42.8\"W"="NA","84º07'25.4\"W"=angle2dec("84º 07' 25.4'' W"),"85º25'05.8 W"=angle2dec("85º 25' 5.8'' W"),"84º07'25.4\"W, 84º08'01.6\"W"="NA",
                                          "84º08'01.6\"W"=angle2dec("84º 08' 01.6'' W"),"84º08'01.6\"W, 85º26'47.33\" W"="NA","84º50'16.7\"W"=angle2dec("84º 50' 16.7'' W")))
unique(data$longitudedecimal)
data$longitudedecimal <- as.numeric(data$longitudedecimal)

#add euthanasia column
data$title <- gsub('[\"]', '', data$title)
unique(sort(data$title))
data$euthanasia<-data$title
data$euthanasia<-revalue(data$euthanasia,
                                c("a coronavirus detected in the vampire bat desmodus rotundus"="Yes, BUT for previous study or rabies surveillance",
                                  "a persistently infecting coronavirus in hibernating myotis lucifugus, the north american little brown bat"="Yes, BUT for previous study or rabies surveillance",
                                  "a real-time pcr assay for bat sars-like coronavirus detection and its application to italian greater horseshoe bat faecal sample surveys"="No",
                                  "alphacoronavirus detected in bats in the united kingdom"="No",
                                  "alphacoronavirus detection in lungs, liver, and intestines of bats from brazil"="Yes, for this study",
                                  "alphacoronavirus in urban molossidae and phyllostomidae bats, brazil"="Yes, BUT for previous study or rabies surveillance",
                                  "alphacoronaviruses detected in french bats are phylogeographically linked to coronaviruses of european bats"="No",
                                  "alphacoronaviruses in new world bats: prevalence, persistence, phylogeny, and potential for interaction with humans"="Yes, BUT for previous study or rabies surveillance",
                                  "amplification of emerging viruses in a bat colony"="No",
                                  "bat coronavirus phylogeography in the western indian ocean"="Yes, BUT for previous study or rabies surveillance",
                                  "bat coronaviruses and experimental infection of bats, the philippines"="Yes, for this study",
                                  "bats are natural reservoirs of sars-like coronaviruses"="No",
                                  "characterization of a new member of alphacoronavirus with unique genomic features in rhinolophus bats"="No",
                                  "circulation of group 2 coronaviruses in a bat species common to urban areas in western europe"="No",
                                  "close relative of human middle east respiratory syndrome coronavirus in bat, south africa"="No",
                                  "coexistence of multiple coronaviruses in several bat colonies in an abandoned mineshaft"="No",
                                  "complete genome sequence of bat coronavirus hku2 from chinese horseshoe bats revealed a much smaller spike gene with a different evolutionary lineage from the rest of the genome"="Yes, for this study",
                                  "coronavirus and paramyxovirus in bats from northwest italy"="No",
                                  "coronavirus antibodies in african bat species"="Yes, for this study",
                                  "coronavirus infection and diversity in bats in the australasian region"="Yes, BUT for previous study or rabies surveillance",
                                  "coronavirus testing indicates transmission risk increases along wildlife supply chains for human consumption in viet nam, 2013-2014"="No",
                                  "coronaviruses detected in bats in close contact with humans in rwanda"="No",
                                  "coronaviruses in bats from mexico"="No",
                                  "coronaviruses in bent-winged bats (miniopterus spp.)"="No",
                                  "coronaviruses in guano from pteropus medius bats in peradeniya, sri lanka"="No",
                                  "coronaviruses in south african bats"="Yes, BUT for previous study or rabies surveillance",
                                  "cross-sectional surveillance of middle east respiratory syndrome coronavirus (mers-cov) in dromedary camels and other mammals in egypt, august 2015 to january 2016"="No",
                                  "detection and characterization of distinct alphacoronaviruses in five different bat species in denmark"="No",
                                  "detection and characterization of diverse alpha- and betacoronaviruses from bats in china"="Yes, for this study",
                                  "detection and full genome characterization of two beta cov viruses related to middle east respiratory syndrome from bats in italy"="No",
                                  "detection and phylogenetic analysis of group 1 coronaviruses in south american bats"="Yes, for this study",
                                  "detection and prevalence patterns of group i coronaviruses in bats, northern germany"="No",
                                  "detection of a virus related to betacoronaviruses in italian greater horseshoe bats"="No",
                                  "detection of alpha and betacoronaviruses in multiple iberian bat species"="No",
                                  "detection of bat coronavirus and specific antibodies in chestnut bat (scotophilus kuhlii) population in central taiwan"="No",
                                  "detection of bat coronaviruses from miniopterus fuliginosus in japan"="Yes, for this study",
                                  "detection of coronavirus genomes in moluccan naked-backed fruit bats in indonesia"="Yes, for this study",
                                  "detection of coronaviruses in bats of various species in italy"="Yes, BUT for previous study or rabies surveillance",
                                  "detection of coronaviruses in pteropus & rousettus species of bats from different states of india"="No, but some bats died during trapping/processing and were necropsied and tissues processed",
                                  "detection of diverse viruses in alimentary specimens of bats in macau"="No",
                                  "detection of group 1 coronaviruses in bats in north america"="No",
                                  "detection of new genetic variants of betacoronaviruses in endemic frugivorous bats of madagascar"="No",
                                  "detection of novel coronaviruses in bats in myanmar"="No",
                                  "detection of novel sars-like and other coronaviruses in bats from kenya"="Yes, for this study",
                                  "detection of polyoma and corona viruses in bats of canada"="Yes, for this study; also rabies submissions + found carcasses were processed",
                                  "detection of potentially novel paramyxovirus and coronavirus viral rna in bats and rats in the mekong delta region of southern viet nam"="No",
                                  "detection of severe acute respiratory syndrome-like, middle east respiratory syndrome-like bat coronaviruses and group h rotavirus in faeces of korean bats"="No",
                                  "detection of specific antibodies to the nucleocapsid protein fragments of severe acute respiratory syndrome-coronavirus and scotophilus bat coronavirus-512 in three insectivorous bat species"="No",
                                  "detection of the severe acute respiratory syndrome-related coronavirus and alphacoronavirus in the bat population of taiwan"="No",
                                  "discovery and genetic analysis of novel coronaviruses in least horseshoe bats in southwestern china"="No",
                                  "discovery of a rich gene pool of bat sars-related coronaviruses provides new insights into the origin of sars coronavirus"="No",
                                  "discovery of novel bat coronaviruses in south china that use the same receptor as middle east respiratory syndrome coronavirus"="No",
                                  "distant relatives of severe acute respiratory syndrome coronavirus and close relatives of human coronavirus 229e in bats, ghana"="No",
                                  "distribution of bat-borne viruses and environment patterns"="Yes, for this study",
                                  "diversity of coronavirus in bats from eastern thailand"="No",
                                  "ecoepidemiology and complete genome comparison of different strains of severe acute respiratory syndrome-related rhinolophus bat coronavirus in china reveal bats as a reservoir for acute, self-limiting infection that allows recombination events"="Yes, for this study",
                                  "evidence for an ancestral association of human coronavirus 229e with bats"="No",
                                  "extensive diversity of coronaviruses in bats from china"="Yes, for this study",
                                  "fatal swine acute diarrhoea syndrome caused by an hku2-related coronavirus of bat origin"="No",
                                  "first report of coronaviruses in northern european bats"="No",
                                  "genetic characteristics of coronaviruses from korean bats in 2016"="No",
                                  "genetic characterization of betacoronavirus lineage c viruses in bats reveals marked sequence divergence in the spike protein of pipistrellus bat coronavirus hku5 in japanese pipistrelle: implications for the origin of the novel middle east respiratory syndrome coronavirus"="Yes, for this study",
                                  "genetic diversity and ecology of coronaviruses hosted by cave-dwelling bats in gabon"="Yes, for this study",
                                  "genetic diversity of bats coronaviruses in the atlantic forest hotspot biome, brazil"="Yes, for this study",
                                  "genetic diversity of coronaviruses in miniopterus fuliginosus bats"="No",
                                  "genomic and serological detection of bat coronavirus from bats in the philippines"="Yes, for this study",
                                  "genomic characterization and infectivity of a novel sars-like coronavirus in chinese bats"="Yes, for this study",
                                  "genomic characterization of severe acute respiratory syndrome-related coronavirus in european bats and classification of coronaviruses based on partial rna-dependent rna polymerase gene sequences"="No",
                                  "group b betacoronavirus in rhinolophid bats, japan"="No",
                                  "human betacoronavirus 2c emc/2012-related viruses in bats, ghana and europe"="No",
                                  "identification of a lineage d betacoronavirus in cave nectar bats (eonycteris spelaea) in singapore and an overview of lineage d reservoir ecology in se asian bats"="No",
                                  "identification of a novel coronavirus in bats"="No",
                                  "identification of alpha and beta coronavirus in wildlife species in france: bats, rodents, rabbits, and hedgehogs"="Yes, BUT for previous study or rabies surveillance",
                                  "identification of diverse alphacoronaviruses and genomic characterization of a novel severe acute respiratory syndrome-like coronavirus from bats in china"="No",
                                  "identification of sars-like coronaviruses in horseshoe bats (rhinolophus hipposideros) in slovenia"="No",
                                  "insectivorous bats carry host specific astroviruses and coronaviruses across different regions in germany"="No",
                                  "intraspecies diversity of sars-like coronaviruses in rhinolophus sinicus and its implications for the origin of sars coronaviruses in humans"="No",
                                  "isolation and characterization of a bat sars-like coronavirus that uses the ace2 receptor"="No",
                                  "longitudinal study of age-specific pattern of coronavirus infection in lyle's flying fox (pteropus lylei) in thailand"="No",
                                  "longitudinal surveillance of betacoronaviruses in fruit bats in yunnan province, china during 2009-2016"="Yes, but only a subset for species ID, tissue tropism work",
                                  "mers-related betacoronavirus in vespertilio superans bats, china"="No",
                                  "middle east respiratory syndrome coronavirus in bats, saudi arabia"="No",
                                  "molecular detection of viruses in kenyan bats and discovery of novel astroviruses, caliciviruses and rotaviruses"="No",
                                  "molecular diversity of coronaviruses in bats"="No",
                                  "molecular identification of betacoronavirus in bats from sardinia (italy): first detection and phylogeny"="No",
                                  "molecular survey of rna viruses in hungarian bats: discovering novel astroviruses, coronaviruses, and caliciviruses"="No",
                                  "neotropical bats from costa rica harbour diverse coronaviruses"="No",
                                  "new alphacoronavirus in mystacina tuberculata bats, new zealand"="No",
                                  "no evidence of coronavirus infection by reverse transcriptase-pcr in bats in belgium"="Yes, BUT for previous study or rabies surveillance",
                                  "novel alphacoronaviruses and paramyxoviruses cocirculate with type 1 and severe acute respiratory system (sars)-related betacoronaviruses in synanthropic bats of luxembourg"="No",
                                  "novel bat alphacoronaviruses in southern china support chinese horseshoe bats as an important reservoir for potential novel coronaviruses"="Yes, for this study",
                                  "novel coronaviruses, astroviruses, adenoviruses and circoviruses in insectivorous bats from northern china"="Yes, for this study",
                                  "persistent infections support maintenance of a coronavirus in a population of australian bats (myotis macropus)"="No",
                                  "prevalence and genetic diversity of coronaviruses in bats from china"="No",
                                  "rapid detection of mers coronavirus-like viruses in bats: pote1ntial for tracking mers coronavirus transmission and animal origin"="N/A (samples from previous study used)",
                                  "recent transmission of a novel alphacoronavirus, bat coronavirus hku10, from leschenault's rousettes to pomona leaf-nosed bats: first evidence of interspecies transmission of coronavirus between bats of different suborders"="Yes, for this study",
                                  "sars-coronavirus ancestor's foot-prints in south-east asian bat colonies and the refuge theory"="No",
                                  "sars-cov related betacoronavirus and diverse alphacoronavirus members found in western old-world"="No",
                                  "seasonal fluctuations of astrovirus, but not coronavirus shedding in bats inhabiting human-modified tropical forests"="No",
                                  "severe acute respiratory syndrome coronavirus-like virus in chinese horseshoe bats"="No",
                                  "surveillance for coronaviruses in bats, lebanon and egypt, 2013-2015"="Yes, but only a subset (no reason provided)",
                                  "surveillance of bat coronaviruses in kenya identifies relatives of human coronaviruses nl63 and 229e and their recombination history"="No",
                                  "the close genetic relationship of lineage d betacoronavirus from nigerian and kenyan straw-colored fruit bats (eidolon helvum) is consistent with the existence of a single epidemiological unit across sub-saharan africa"="No",
                                  "the persistent prevalence and evolution of cross-family recombinant coronavirus gccdc1 among a bat population: a two-year follow-up"="No",
                                  "viral diversity of bat communities in human-dominated landscapes in mexico"="No",
                                  "viral diversity of microbats within the south west botanical province of western australia"="No",
                                  "virome analysis for identification of novel mammalian viruses in bats from southeast china"="Yes, for this study",
                                  "virus survey in populations of two subspecies of bent-winged bats (miniopterus orianae bassanii and oceanensis) in south-eastern australia reveals a high prevalence of diverse herpesviruses"="No, but some bats were found dead and processed for study",
                                  "white-nose syndrome is associated with increased replication of a naturally persisting coronaviruses in bats"="N/A; experimental infection",
                                  "wide diversity of coronaviruses in frugivorous and insectivorous bat species: a pilot study in guinea, west africa"="No",
                                  "wild animal surveillance for coronavirus hku1 and potential variants of other coronaviruses"="No"))

#simplify study detection method (specific)
unique(data$method_specific)
data$method_specific[which(data$method_specific=="")] <- "NA"
data$method_specific_simplified <- data$method_specific
data$method_specific_simplified <- revalue(data$method_specific_simplified, c("Lateral flow immunoassay"="Lateral Flow Immunoassay","broadly reactive heminested RT-PCR"="PCR","hemi-nested RT-PCR"="PCR","RT-PCR"="PCR","RT-PCR,nested RT-PCR"="PCR","semi-nested RT-PCR"="PCR","nested RT-PCR"="PCR","nested PCR"="PCR","end point PCR"="PCR","RT-PCR,hemi-nested PCR"="PCR","consensus RT-PCR"="PCR","consensus nested-PCR"="PCR","\"SARS coronavirus crude antigen ELISA (Yu, 2008)\""="ELISA","double heminested RT-PCR"="PCR","rtRT-PCR"="PCR","RT-PCR,semi-nested RT-PCR,nested RT-PCR"="PCR","broadly reactive nested RT-PCR"="PCR","pancoronavirus Nested PCR"="PCR","pancoronavirus nested RT-PCR"="PCR","degenerate concensus PCR"="PCR","nested RT-PCR,semi-nested RT-PCR"="PCR","Enzyme immunoassay"="ELISA","RT-qPCR,semi-nested PCR"="PCR","pancoronavirus RT-PCR"="PCR"))
unique(data$method_specific_simplified)

#check virus genera data
unique(data$virus_genus)
data$virus_genus[which(data$virus_genus=="betacoronavirus ")] <- "betacoronavirus"

##load in Upham phylogeny
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

## are all bats in tree
setdiff(data$bat_species,tree$tip)

## multiple species?
data$spec=sapply(strsplit(data$bat_species,","),function(x) length(x))
which(data$spec==2)
data$bat_species[which(data$spec==2)]
#only "Peking Myotis, Myotis pequinius" + pooled rhinopholus + myotis species - fixed below 

#fix bats for dataset being given to reader (e.g. rows with only genus-level info are kept)
data$species_for_reader=data$bat_species
#species not in phylogeny, not enough info (e.g. Insectivorous bats) -> Drop 
data$species_for_reader=revalue(data$species_for_reader,
                     c("Artibeus phaeotis"="Dermanura phaeotis",
                       "Artibeus literatus"="Artibeus lituratus",
                       "Vampyressa pussilla"="Vampyressa pusilla",
                       "Myotis cillolabrum"="Myotis ciliolabrum",
                       "Macronycteris commersoni"="Hipposideros commersoni",
                       "Chaerephon leucogaster"="Chaerephon pumilus",
                       "Momopterus acetabulosus"="Mormopterus acetabulosus",
                       "Taphosius mauritianus"="Taphozous mauritianus",
                       "Chaerephon sp"="Chaerephon",
                       "Chaerephon pusillus"="Chaerephon pumilus",
                       "Pteropus seychellensis comorensis"="Pteropus seychellensis",
                       "Rhinolophus lobatus"="Rhinolophus landeri",
                       "Rhinolophus rhodesiae"="Rhinolophus simulator",
                       "Rhinolophus sp"="Rhinolophus",
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
                       "Rhinolophus spp"="Rhinolophus",
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
                       "Artibeus watsoni"="Dermanura watsoni",
                       "Mormoops megalohyla"="Mormoops megalophylla",
                       "Pteropus medius"="Pteropus giganteus",
                       "Myonicteris angolensis"="Myonycteris angolensis",
                       "Stenonycteris lanosus"="Rousettus lanosus",
                       "Pteropus sp"="Pteropus",
                       "Microchiroptera (order: chiroptera)"="Microchiroptera",
                       "Hipposideros terasensis"="Hipposideros armiger",
                       "Scotorepens rueppellii"="Scoteanax rueppellii",
                       "Cynopterus spp"="Cynopterus",
                       "Eonycteris spp"="Eonycteris",
                       "Macroglossus spp"="Macroglossus",
                       "Rousettus spp"="Rousettus",
                       "Taphozous spp"="Taphozous",
                       "Chalinolobus spp"="Chalinolobus",
                       "Scotophilus spp"="Scotophilus",
                       "Scotorepens spp"="Scotorepens",
                       "Chaerephon spp"="Chaerephon",
                       "Molossidae bats"="Mollosidae",
                       "Neuromicia nana"="Neoromicia nana",
                       "Neuromicia helios"="Neoromicia helios",
                       "Neoromicia spp"="Neoromicia",
                       "Nycticeinops schlieffenii"="Nycticeinops schlieffeni",
                       "Rhinolophus darlingi damarensis"="Rhinolophus darlingi",
                       "Insectivorous bat species"="Drop",
                       "Rousettous aegyptiacus"="Rousettus aegyptiacus",
                       "Taphozous perforates"="Taphozous perforatus",
                       "Coelops frithii formosanus"="Coelops frithii",
                       "Hipposideros armiger terasensis"="Hipposideros armiger",
                       "Barbastella darjelingesis"="Barbastella beijingensis",
                       "Myotis fimbriatus taiwanensis"="Myotis fimbriatus",
                       "Myotis rufoniger"="Drop",
                       "Pipistrellus montanus"="Drop",
                       "Pipistrellus taiwanesis"="Drop",
                       "Myotis formosus flavus"="Myotis formosus",
                       "Fruit bat"="Drop",
                       "Chaerephon spp."="Chaerephon",
                       "Pipistrellus nanus"="Neoromicia nana",
                       "Hipposideros caffer ruber"="Hipposideros ruber",
                       "Myotis sp."="Myotis",
                       "Rhinolophus rouxi"="Rhinolophus rouxii",
                       "Pipistrellus sp."="Pipistrellus",
                       "Rousettus leschenaultia"="Rousettus leschenaultii",
                       "Vespertilio superans"="Vespertilio sinensis",
                       "Rhinolophus affiinus"="Rhinolophus affinis",
                       "Rousettus lechenaulti"="Rousettus leschenaultii",
                       "Tyloncyteris pachypus"="Tylonycteris pachypus",
                       "Hypsugo pulveratus"="Pipistrellus pulveratus",
                       "Miniopterus filiginosus"="Miniopterus fuliginosus",
                       "Peking Myotis,Myotis pequinius"="Myotis pequinius",
                       "Rhinolophus malyanus"="Rhinolophus malayanus",
                       "Coelops frithi"="Coelops frithii",
                       "Craseonycteris thonlongyal"="Craseonycteris thonglongyai",
                       "Rhinoloplus malayanus"="Rhinolophus malayanus",
                       "Rhinoloplus sp."="Rhinolophus",
                       "Miniopterus meghrebensis"="Miniopterus maghrebensis",
                       "Rhinolophus borneenis"="Rhinolophus borneensis",
                       "Tadarida spp"="Tadaria",
                       "Myotis spp"="Myotis",
                       "Pipistrellus minus"="Pipistrellus tenuis",
                       "Pipistrellus spp"="Pipistrellus",
                       "Tylonycteris spp"="Tylonycteris",
                       "Hipposideros cf. gigas"="Hipposideros gigas",
                       "Hipposideros cf. ruber"="Hipposideros ruber",
                       "Nycteris cf. gambiensis"="Nycteris gambiensis",
                       "Lissonycteris angolensis"="Myonycteris angolensis",
                       "Murina spp"="Murina",
                       "Megaerops kusnotei"="Megaerops kusnotoi",
                       "Rhinolophus affinus"="Rhinolophus affinis",
                       "Molossus major"="Molossus molossus",
                       "Mormoops sp."="Mormoops",
                       "Dobsonia mollucensis"="Dobsonia moluccensis",
                       "Haplonicteris fischeri"="Haplonycteris fischeri",
                       "Miniopterus orianae bassanii"="Miniopterus schreibersii",
                       "Miniopterus orianae oceanensis"="Miniopterus schreibersii",
                       "Nyctiphilus geoffroyi"="Nyctophilus geoffroyi",
                       "Nyctophilus major"="Nyctophilus timoriensis",
                       "Austronomus australis"="Tadarida australis",
                       "Ozimops sp"="Ozimops",
                       "Hipposideros sp"="Hipposideros",
                       "Nycteris sp."="Nycteris",
                       "Epomps buettikoferi"="Epomops buettikoferi",
                       "Myotis oxygnatus"="Myotis blythii",
                       "Myotis capaccini"="Myotis capaccinii",
                       "Rhinolophus cornutus"="Drop",
                       "Pipistrelus abramus"="Pipistrellus abramus",
                       "Myotis formosus chofukusei"="Myotis formosus",
                       "Pteropus spp."="Pteropus",
                       "Taphozous"="Taphozous",
                       "Miniopterus africanus"="Miniopterus inflatus",
                       "Vampyriscus nymphaea"="Vampyressa nymphaea",
                       "Glossophaga comissarisi"="Glossophaga commissarisi",
                       "Eptesicus furnalis"="Eptesicus furinalis",
                       "Miniopterus sp."="Miniopterus",
                       "Scotoecus sp."="Miniopterus",
                       "Rhinolophus sinicus,Rhinolophus ferrumequinum"="Rhinolophus",
                       "Myotis aurascens Kuzyakin,Myotis petax"="Myotis",
                       "(multiple) species of horseshoe bat"="Rhinolophus"))

## check
unique(data$species_for_reader)
data_for_reader <- data

#remove genus-only/drop rows for analyses (trim dataset to only phylogeny) 
setdiff(data$species_for_reader,tree$tip)
set=data[data$species_for_reader%in%tree$tip.label,,]
set$species=factor(set$species_for_reader)
set$studies=factor(set$title)

###genus-level coronavirus datasets: 
#1) infection prevalence analyses (antibody + IHC rows removed, flagged rows removed)
set_infection_prev <- set[which(set$detection_method!= "IHC" | set$detection_method!="antibody"),,]
set_infection_prev <- set_infection_prev[which(set_infection_prev$Flag=="only use these two rows for prevalence estimates" | set_infection_prev$Flag==""),,]
#^dataset for prevalence estimates

#2) other analyses (two rows only for prevalence have been removed)
set_other <- set[which(set$Flag=="" | set$Flag=="include rows for taxonomic/geographic patterns in yes/no sampled and yes/no positive but NOT in prevalence analyses because denominators correspond to region and not species" | set$Flag=="airtable double counts by including study table 1 (ID 1120-1122)+ table 2 (ID 2945-2956): study included in binary analyses but not prevalence estimates"),,]

###alpha-only and beta-only datasets:
#1) infection prevalence analyses
unique(set_infection_prev$virus_genus)
set_infection_prev_alphaonly <- set_infection_prev[which(set_infection_prev$virus_genus=="alphacoronavirus" | set_infection_prev$virus_genus=="alphacoronavirus/alphacoronavirus coinfection"),,] 
set_infection_prev_betaonly <- set_infection_prev[which(set_infection_prev$virus_genus=="betacoronavirus"),,]

#2) other analyses
set_other_alphaonly <- set_other[which(set_other$virus_genus=="alphacoronavirus" | set_other$virus_genus=="alphacoronavirus/alphacoronavirus coinfection"),,] 
set_other_betaonly <- set_other[which(set_other$virus_genus=="betacoronavirus"),,]

#infection prevalence all coronaviruses
## trim tree to species in set
stree=keep.tip(tree,as.character(unique(set_infection_prev$species)))

## convert tree to correlation matrix
cmatrix=vcv.phylo(stree,cor=T)

## make observation and study-level random effect
set_infection_prev$observation=factor(1:nrow(set_infection_prev))
set_infection_prev$study=factor(set_infection_prev$title)

## pft in escalc for yi and vi 
set_infection_prev=data.frame(set_infection_prev,escalc(xi=set_infection_prev$positives,ni=set_infection_prev$sample,measure="PFT"))

## back transform
set_infection_prev$backtrans=transf.ipft(set_infection_prev$yi,set_infection_prev$sample)

## check
plot(set_infection_prev$prevalence,set_infection_prev$backtrans)
abline(0,1)

## species and phylo effect
set_infection_prev$phylo=set_infection_prev$species
set_infection_prev$species=set_infection_prev$phylo

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
set_infection_prev$species_char <- as.character(set_infection_prev$species)

which(taxa$species==set_infection_prev$species_char[1])

data_fam <- data.frame(family = taxa$fam[986])

nrow(set_infection_prev)
for(i in 2:2281) {
  new_row <- data.frame(family = taxa$fam[which(taxa$species==set_infection_prev$species_char[i])])
  data_fam <- rbind(data_fam, new_row)
}

set_infection_prev <- cbind(set_infection_prev, data_fam)



#export to github desktop folder
#reader
write.csv(data_for_reader,"~/Documents/GitHub/batgap/data/data_for_reader.csv")
#infection prevalence
write.csv(set_infection_prev,"~/Documents/GitHub/batgap/data/set_infection_prevalence.csv")
write.csv(set_infection_prev_alphaonly,"~/Documents/GitHub/batgap/data/set_infection_prevalence_alphaonly.csv")
write.csv(set_infection_prev_betaonly,"~/Documents/GitHub/batgap/data/set_infection_prevalence_betaonly.csv")
#other analyses
write.csv(set_other,"~/Documents/GitHub/batgap/data/set_other.csv")
write.csv(set_other_alphaonly,"~/Documents/GitHub/batgap/data/set_other_alphaonly.csv")
write.csv(set_other_betaonly,"~/Documents/GitHub/batgap/data/set_other_betaonly.csv")



#ready for models





