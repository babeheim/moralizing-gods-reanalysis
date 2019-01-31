#setwd("/Users/pesavage/Documents/Research/Oxford Seshat/Data/MoralizingGods")

map<-read.csv("map.csv",header=TRUE,row.names=1)
x<-map$Longitude
y<-map$Latitude
z<-2.8*sqrt(map$SC)

#USING MAPS


#all NGAs outlined coloured by moralising gods, area proportional to SC, text = kya

map$kya[20]=".2"
map$kya[19]=".2"
map$kya[18]=".3"
map$kya[17]=".5"
map$kya[16]=".9"

map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0)) 
points(subset(x,map$Religion=="Zoroastrianism"),subset(y,map$Religion =="Zoroastrianism"), bg="lightcoral", pch=21,cex=subset(z,map$Religion =="Zoroastrianism"))
points(subset(x,map$Religion=="Abrahamic"),subset(y,map$Religion =="Abrahamic"), bg="orange", pch=21,cex=subset(z,map$Religion =="Abrahamic"))
points(subset(x,map$Religion=="Buddhism"),subset(y,map$Religion =="Buddhism"), bg="dodgerblue", pch=21,cex=subset(z,map$Religion =="Buddhism"))
points(subset(x,map$Religion=="OtherMHG"),subset(y,map$Religion =="OtherMHG"), bg="yellow", pch=21,cex=subset(z,map$Religion =="OtherMHG"))
points(subset(x,map$Religion =="OtherBSP"),subset(y,map$Religion =="OtherBSP"), bg="mediumpurple1", pch=21,cex=subset(z,map$Religion =="OtherBSP"))
points(subset(x,map$Religion =="N"),subset(y,map$Religion =="N"), bg="grey", pch=21,cex=subset(z,map$Religion =="N"))

text(x,y, labels = map$kya, cex=.55,col="black")

pos<-legend(-182,1, c("Zoroastrianism","Abrahamic","Other MHG","Buddhism","Other BSP","Absent"),pch=21, pt.bg=c('lightcoral', 'orange',"yellow","dodgerblue","mediumpurple1","grey"), bty='n', cex=0.75,pt.cex=1.6,title=expression(paste(bold("Earliest precolonial evidence\nof moralizing gods (kya)"))))
xleft <- pos$rect[["left"]]
ytop <- pos$rect[["top"]]
ybottom <- ytop - pos$rect[["h"]]
xright <- xleft + pos$rect[["w"]]
rect(xleft, ybottom+2, xright, ytop+8)
