# Pat Savage (2017, University of Oxford)

#setwd("/Users/pesavage/Documents/Research/Oxford Seshat/Data/SCBigGodsOct2017")

#library("ggmap") 
#library(maptools) 
#library(maps) 

map<-read.csv("map.csv",header=TRUE,row.names=1)
x<-map$Longitude
y<-map$Latitude
z<-3*sqrt(map$SC)

cex<-1.9

#USING MAPS

#all NGAs outlined coloured by moralising gods, area proportional to SC, no distinctino between BG and BSP

map("world", fill=TRUE, col="white", bg="gray90", ylim=c(-60, 90), mar=c(0,0,0,0)) 
points(subset(x,map$MG=="Y"),subset(y,map$MG=="Y"), bg="gray50", pch=21,cex=subset(z,map$MG=="Y"))
points(subset(x,map$MG=="N"),subset(y,map$MG=="N"), bg="gray90", pch=21,cex=subset(z,map$MG=="N"))

text(x,y, labels = row.names(map), cex=.6,col="black")

legend(-180,-23, c("Present","Absent"),pch=21, pt.bg=c('gray50',"gray90"), bty='y', cex=0.7,pt.cex=1.6,title="Precolonial evidence of moralizing gods")