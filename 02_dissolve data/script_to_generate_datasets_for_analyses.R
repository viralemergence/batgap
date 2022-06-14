#again, need tree + phylogeny
#download tree and Upham phylogeny from github.com/viralemergence/batgap/01-generate_datacov
#tree: "MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre"
tree <- read.nexus("MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre")
#phylogeny: "taxonomy_mamPhy_5911species.csv"
taxa <- read.csv('taxonomy_mamPhy_5911species.csv',header=T)
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

#download datacov from https://github.com/viralemergence/datacov
datacov <- read.csv("datacov.csv")

#remove genus-only/drop rows for analyses (trim dataset to only phylogeny) 
setdiff(datacov$species_for_reader,tree$tip)
dataset=datacov[datacov$species_for_reader%in%tree$tip.label,,]
dataset$species=factor(dataset$species_for_reader)
dataset$studies=factor(dataset$title)

###datasets: 
#1) infection prevalence analyses (antibody + IHC rows removed, flagged rows removed, coinfection rows removed) - all cov genera
set_infection_prev <- dataset[-which(dataset$detection_method== "IHC" | dataset$detection_method=="antibody"),,]
set_infection_prev <- set_infection_prev[which(set_infection_prev$Flag=="only use these two rows for prevalence estimates" | set_infection_prev$Flag==""),,]
set_infection_prev <- set_infection_prev[-which(set_infection_prev$virus_genus_simplified== "coinfection"),,]

#2) other analyses (five rows only for prevalence have been removed, coinfection rows removed) - all cov genera
set_other <- dataset[which(dataset$Flag=="" | dataset$Flag=="include rows for taxonomic/geographic patterns in yes/no sampled and yes/no positive but NOT in prevalence analyses because denominators correspond to region and not species" | dataset$Flag=="airtable double counts by including study table 1 (ID 1120-1122)+ table 2 (ID 2945-2956): study included in binary analyses but not prevalence estimates"),,]
set_other <- set_other[-which(set_other$virus_genus_simplified== "coinfection"),,]

###alpha-only and beta-only datasets:
#1) infection prevalence analyses
set_infection_prev_alphaonly <- set_infection_prev[c(which(set_infection_prev$virus_genus_simplified=="alphacoronavirus"),which(set_infection_prev$virus_genus_simplified=="genus not specified" & set_infection_prev$positives==0)),,] 
set_infection_prev_betaonly <- set_infection_prev[c(which(set_infection_prev$virus_genus_simplified=="betacoronavirus"),which(set_infection_prev$virus_genus_simplified=="genus not specified" & set_infection_prev$positives==0)),,]

#2) other analyses
set_other_alphaonly <- set_other[c(which(set_other$virus_genus_simplified=="alphacoronavirus"),which(set_other$virus_genus_simplified=="genus not specified" & set_other$positives==0)),,] 
set_other_betaonly <- set_other[c(which(set_other$virus_genus_simplified=="betacoronavirus"),which(set_other$virus_genus_simplified=="genus not specified" & set_other$positives==0)),,]

##add column for bat family 
#all genera
set_infection_prev$species_char <- as.character(set_infection_prev$species)
which(taxa$species==set_infection_prev$species_char[1])
data_fam <- data.frame(family = taxa$fam[986])
nrow(set_infection_prev)
for(i in 2:2076) {
  new_row <- data.frame(family = taxa$fam[which(taxa$species==set_infection_prev$species_char[i])])
  data_fam <- rbind(data_fam, new_row)
}
set_infection_prev <- cbind(set_infection_prev, data_fam)
#alpha only
set_infection_prev_alphaonly$species_char <- as.character(set_infection_prev_alphaonly$species)
which(taxa$species==set_infection_prev_alphaonly$species_char[1])
data_fam <- data.frame(family = taxa$fam[986])
nrow(set_infection_prev_alphaonly)
for(i in 2:1770) {
  new_row <- data.frame(family = taxa$fam[which(taxa$species==set_infection_prev_alphaonly$species_char[i])])
  data_fam <- rbind(data_fam, new_row)
}
set_infection_prev_alphaonly <- cbind(set_infection_prev_alphaonly, data_fam)
#beta only
set_infection_prev_betaonly$species_char <- as.character(set_infection_prev_betaonly$species)
which(taxa$species==set_infection_prev_betaonly$species_char[1])
data_fam <- data.frame(family = taxa$fam[986])
nrow(set_infection_prev_betaonly)
for(i in 2:1627) {
  new_row <- data.frame(family = taxa$fam[which(taxa$species==set_infection_prev_betaonly$species_char[i])])
  data_fam <- rbind(data_fam, new_row)
}
set_infection_prev_betaonly <- cbind(set_infection_prev_betaonly, data_fam)

#infection prevalence datasets
write.csv(set_infection_prev,"~/Documents/GitHub/batgap/02a-datasets_for_analyses/set_infection_prevalence.csv")
write.csv(set_infection_prev_alphaonly,"~/Documents/GitHub/batgap/02a-datasets_for_analyses/set_infection_prevalence_alphaonly.csv")
write.csv(set_infection_prev_betaonly,"~/Documents/GitHub/batgap/02a-datasets_for_analyses/set_infection_prevalence_betaonly.csv")

#other analyses datasets
write.csv(set_other,"~/Documents/GitHub/batgap/02a-datasets_for_analyses/set_other.csv")
write.csv(set_other_alphaonly,"~/Documents/GitHub/batgap/02a-datasets_for_analyses/set_other_alphaonly.csv")
write.csv(set_other_betaonly,"~/Documents/GitHub/batgap/02a-datasets_for_analyses/set_other_betaonly.csv")