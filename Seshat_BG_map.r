# Pat Savage (2017, University of Oxford)

#setwd("/Users/pesavage/Documents/Research/Oxford Seshat/Data/SCBigGodsOct2017")

library("ggmap") 
library(maptools) 
library(maps) 

map<-read.csv("map.csv",header=TRUE,row.names=1)
x<-map$Longitude
y<-map$Latitude
z<-3*sqrt(map$SC)

cex<-1.9
#USING MAPS
#points only 
#map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0)) 
#points(x,y, bg="red", pch=21)

#grayscale 
#map("world", fill=TRUE, col="white", bg="gray90", ylim=c(-60, 90), mar=c(0,0,0,0)) 
#points(x,y, bg=c(rep("gray33",19),rep("gray66",16)), pch=21,cex=1.5)

#USING MAPS, labeled/numbered
#only 12 pre/post NGAs outlined
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0)) 
points(subset(x,map$SupernaturalPunishment=="MHG"),subset(y,map$SupernaturalPunishment=="MHG"), col="red", pch=19,cex=cex)
points(subset(x,map$SupernaturalPunishment=="BSP"),subset(y,map$SupernaturalPunishment=="BSP"), col="orange", pch=19,cex=cex)
points(subset(x,map$SupernaturalPunishment=="N"),subset(y,map$SupernaturalPunishment=="N"), col="grey", pch=19,cex=cex)
points(subset(x,map$PrePostSPData=="Y"),subset(y, map$PrePostSPData=="Y"), col="black", pch=1,cex=cex,lwd=2)

text(x,y, labels = row.names(map), cex=.6,col="black")

legend(-180,-25, c("Big Gods","Broad Supernatural Punishment","None"),pch=21, pt.bg=c('red', 'orange',"grey"), bty='y', cex=0.7,pt.cex=1.6,title="Earliest evidence of moralising gods")
legend(-45,-32, c("Before and after appearance of moralising gods","Pre-/post-moralising gods only'"),pch=21, pt.bg="grey", lwd=c(3,0), bty='y', cex=0.7,pt.cex=1.6,title="Social complexity data availability")

#all NGAs outlined gray only
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0)) 
points(x,y, bg="grey", pch=21,cex=cex)

#text(x,y, labels = row.names(map), cex=.6,col="black")

#all NGAs outlined coloured by moralising gods
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0)) 
points(subset(x,map$SupernaturalPunishment=="MHG"),subset(y,map$SupernaturalPunishment=="MHG"), bg="red", pch=21,cex=cex)
points(subset(x,map$SupernaturalPunishment=="BSP"),subset(y,map$SupernaturalPunishment=="BSP"), bg="orange", pch=21,cex=cex)
points(subset(x,map$SupernaturalPunishment=="N"),subset(y,map$SupernaturalPunishment=="N"), bg="grey", pch=21,cex=cex)

text(x,y, labels = row.names(map), cex=.6,col="black")

legend(-180,-25, c("Big Gods","Broad Supernatural Punishment","None"),pch=21, pt.bg=c('red', 'orange',"grey"), bty='y', cex=0.7,pt.cex=1.6,title="Earliest evidence of moralizing gods")

#all NGAs outlined coloured by moralising gods, area proportional to SC

map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0)) 
points(subset(x,map$SupernaturalPunishment=="MHG"),subset(y,map$SupernaturalPunishment=="MHG"), bg="red", pch=21,cex=subset(z,map$SupernaturalPunishment=="MHG"))
points(subset(x,map$SupernaturalPunishment=="BSP"),subset(y,map$SupernaturalPunishment=="BSP"), bg="orange", pch=21,cex=subset(z,map$SupernaturalPunishment=="BSP"))
points(subset(x,map$SupernaturalPunishment=="N"),subset(y,map$SupernaturalPunishment=="N"), bg="grey", pch=21,cex=subset(z,map$SupernaturalPunishment=="N"))

text(x,y, labels = row.names(map), cex=.6,col="black")

legend(-180,-25, c("Big Gods","Broad Supernatural Punishment","None"),pch=21, pt.bg=c('red', 'orange',"grey"), bty='y', cex=0.7,pt.cex=1.6,title="Earliest evidence of moralizing gods")

#all NGAs outlined coloured by moralising gods, area proportional to SC, no distinctino between BG and BSP

map("world", fill=TRUE, col="white", bg="gray90", ylim=c(-60, 90), mar=c(0,0,0,0)) 
points(subset(x,map$MG=="Y"),subset(y,map$MG=="Y"), bg="gray50", pch=21,cex=subset(z,map$MG=="Y"))
points(subset(x,map$MG=="N"),subset(y,map$MG=="N"), bg="gray90", pch=21,cex=subset(z,map$MG=="N"))

text(x,y, labels = row.names(map), cex=.6,col="black")

legend(-180,-23, c("Present","Absent"),pch=21, pt.bg=c('gray50',"gray90"), bty='y', cex=0.7,pt.cex=1.6,title="Precolonial evidence of moralizing gods")
