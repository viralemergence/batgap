## bat coronavirus gap analysis
## 01_lit search clean
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## libraries
library(Hmisc)
library(gtools)
library(stringr)
library(XML)
library(plyr)

## load searches
setwd("~/Desktop/batgap/raw searches")

## wos
wos=read.csv("savedrecs wos.csv")

## trim
keeps=c("Title","Authors","Source.Title","Publication.Year","DOI")
wos=wos[keeps]; rm(keeps)

## clean data
wdata=wos

## fix names
names(wdata)=c("title","author","journal","year","doi")

## save old author
wdata$author2=wdata$author

## assign study ID
wdata$ID=factor(1:nrow(wdata),levels=1:nrow(wdata))

## function to fix name
nfix=function(x){
  
  ## subset
  set=wdata[which(wdata$ID==x),]
  
  ## split author list
  auth=strsplit(set$author,"; ")
  
  ## split into last names
  auth=strsplit(auth[[1]],",")
  auth=sapply(auth,function(x) x[1])
  
  ## lower case everything, then upper case first
  auth=tolower(auth)
  auth=capitalize(auth)
  
  ## get only last names
  auth=paste(as.character(auth),collapse="; ")
  
  ## split again
  auth2=strsplit(auth,"; ")
  
  ## get length
  len=length(auth2[[1]])
  
  ## if else
  name=ifelse(len==1,auth2[[1]],
              ifelse(len==2,paste(auth2[[1]][1],"&",auth2[[1]][2]),paste(auth2[[1]][1],"et al.")))
  return(name)
}

## fix name
wdata$author=sapply(levels(wdata$ID),nfix)

## clean
wdata$author2=NULL

## capitalize first letter
simpleCap=function(x){
  
  ## lower case everything
  s=tolower(x)
  
  ## split apart
  s=strsplit(s," ")[[1]]
  
  ## paste
  s=ifelse(s%in%c("a","of","to","and"),s,capitalize(s))
  s=paste(s,collapse=" ")
  return(s)
}

## fix journal
wdata$journal=sapply(wdata$journal,simpleCap)

## make study
wdata$study=with(wdata,paste(author,year,journal))

## add source
wdata$source="WOS"

## fix title to lower
wdata$title=tolower(wdata$title)

## clean
rm(wos)

## pubmed
pubmed=read.csv("csv-batORChiro-set.csv")

## trim
pdata=pubmed[c("Title",'Authors',"Journal.Book","Publication.Year",'DOI')]

## fix names
names(pdata)=c("title","author","journal","year","doi")

## save old author
pdata$author2=pdata$author

## assign study ID
pdata$ID=factor(1:nrow(pdata),levels=1:nrow(pdata))

## function to fix name
nfix=function(x){
  
  ## subset
  set=pdata[which(pdata$ID==x),]
  
  ## split author list
  auth=strsplit(set$author,", ")
  
  ## split into last names
  auth=strsplit(auth[[1]]," ")
  auth=sapply(auth,function(x) x[1])
  
  ## lower case everything, then upper case first
  auth=tolower(auth)
  auth=capitalize(auth)
  
  ## get only last names
  auth=paste(as.character(auth),collapse="; ")
  
  ## split again
  auth2=strsplit(auth,"; ")
  
  ## get length
  len=length(auth2[[1]])
  
  ## if else
  name=ifelse(len==1,auth2[[1]],
              ifelse(len==2,paste(auth2[[1]][1],"&",auth2[[1]][2]),paste(auth2[[1]][1],"et al.")))
  return(name)
}

## fix name
pdata$author=sapply(levels(pdata$ID),nfix)

## clean
pdata$author2=NULL

## fix journal
pdata$journal=sapply(pdata$journal,simpleCap)

## make study
pdata$study=with(pdata,paste(author,year,journal))

## add source
pdata$source="PubMed"

## fix title to lower
pdata$title=tolower(pdata$title)

## clean
rm(pubmed)

## load global health
cdata=read.csv("cab clean.csv")

## trim 
cdata=cdata[c("journal","year","auth1","auth2","auth3","auth4","title","doi")]
cdata$source="CAB"

## add id
cdata$ID=factor(1:nrow(cdata),levels=1:nrow(cdata))

## fix to lower
cdata$title=tolower(cdata$title)

## fix title
premove=function(x){
  
  ## subset
  set=cdata[which(cdata$ID==x),]
  
  ## get title
  tit=as.character(set$title)
  
  ## separate all
  tit=paste(strsplit(tit,"[.]")[[1]],collapse="")
  return(tit)
}
cdata$title=sapply(levels(cdata$ID),premove)

## combine author columns
cdata$author=with(cdata,paste(auth1,auth2,auth3,auth4,sep="; "))

## clean old columns
cdata$auth1=NULL
cdata$auth2=NULL
cdata$auth3=NULL
cdata$auth4=NULL

## function to fix name
nfix=function(x){
  
  ## subset
  set=cdata[which(cdata$ID==x),]
  
  ## split author list
  auth=strsplit(set$author,"; ")
  
  ## split into last names
  auth=strsplit(auth[[1]],",")
  auth=sapply(auth,function(x) x[1])
  auth=auth[!is.na(auth)]
  
  ## get length
  len=length(auth)
  
  ## if else
  name=ifelse(len==1,auth,
              ifelse(len==2,paste(auth[1],"&",auth[2]),paste(auth[1],"et al.")))
  return(name)
}

## fix authors
cdata$author2=sapply(levels(cdata$ID),nfix)
cdata$author=cdata$author2
cdata$author2=NULL

## make study variable
cdata$study=with(cdata,paste(author,year,journal,sep=" "))

## clean function
rm(nfix,simpleCap,premove)

## match titles
cols=names(wdata)[names(wdata)%in%names(pdata)]
wdata=wdata[cols]
pdata=pdata[cols]
cdata=cdata[cols]
rm(cols)

## combine
lit_data=smartbind(wdata,pdata,cdata)

## clean
rm(wdata,pdata,cdata)

## clean titles
lit_data$title=str_replace(lit_data$title,'\\[','')
lit_data$title=str_replace(lit_data$title,'\\]','')

## sort
lit_data=lit_data[order(lit_data$title),]

## remove
lit_data2=lit_data[!duplicated(lit_data$title),]

## save good doi
lit_data2$doi_pure=lit_data2$doi

## make dummy doi
lit_data2$doi2=1:nrow(lit_data2)
lit_data2$doi=ifelse(lit_data2$doi=='',lit_data2$doi2,lit_data2$doi)
lit_data2$doi2=NULL

## ditto for doi
lit_data2=lit_data2[!duplicated(lit_data2$doi),]

## fix doi
lit_data2$doi=lit_data2$doi_pure
lit_data2$doi_pure=NULL

## tally
nrow(lit_data)-nrow(lit_data2)

## clean
lit_data=lit_data2
rm(lit_data2)

## new id
lit_data$ID=1:nrow(lit_data)

## add columns
lit_data$exclude_title=""
lit_data$exclude_reason=''
lit_data$include=''
lit_data$entered=''

## export
setwd("~/Desktop/batgap/data")
write.csv(lit_data,'cleaned initial batgap search.csv')