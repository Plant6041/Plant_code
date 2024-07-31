library(sp)
library(raster)
library(foreach)
library(iterators)
library(parallel)
library(doParallel)
library(gdm)
library(dplyr)

##import species and environmental data
sppData <- read.csv("example_species_data.csv")
envTab <- read.csv("example_environment_data.csv")

##format site pair using species and environmental data
sitePairTab <- formatsitepair(sppData, 2, XColumn="Long", YColumn="Lat", sppColumn="species", 
                              siteColumn="site", predData=envTab)

##fit GDM
gdmTabMod <- gdm(sitePairTab, geo=TRUE)
summary(gdmTabMod)
plot(gdmTabMod, plot.layout=c(2,4))

##import habitat condition data
hj2018 <- read.csv("example_hj2018.csv")
hj2018 <- hj2018$hj2018

##reorganize the data to run the SAR loop
defaultdistance <- c(1)
defaultweight <- c(1)
n=5
rtpenvrast <- envTab[c(2:8)]
sitepi <- rtpenvrast[0,1:6]
colnames(sitepi) <- c("s1.xCoord","s1.yCoord","sumsij","sumhj2018sij","persistence","extinction")

##run the SAR loop
for (i in 1:nrow(rtpenvrast)){
  sitei <- rtpenvrast[i,]
  rtpenvrasttemp <- (rtpenvrast[-i,])
  hj2018temp <- (hj2018[-i])
  sitepair <- cbind(rtpenvrasttemp,sitei,defaultdistance,defaultweight)[c((2*n+5),(2*n+6),(n+3),(n+4),1,2,(n+5):(2*n+4),3:(n+2))]
  colnames(sitepair) <- colnames(sitePairTab)
  dij <- predict(gdmTabMod,sitepair)
  sij <- (1-dij)
  sitepairisij <- cbind(sitepair[,c(3:6)],sij,hj2018temp)
  
  sitegroup <- sitepairisij%>%
    group_by(s1.xCoord,s1.yCoord)%>%
    summarise(sumsij=sum(sij),
              sumhj2018sij=sum(hj2018temp*sij))
  sitegroup$persistence <- (sitegroup$sumhj2018sij/sitegroup$sumsij)^0.25
  sitegroup$extinction <- 1-sitegroup$persistence
  sitepi <- rbind(sitepi,sitegroup[,c(1:6)])
}

##export result
write.csv(sitepi,file = "example_extinction_risk.csv")
