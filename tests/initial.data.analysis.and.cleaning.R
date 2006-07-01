# source("CS.mpp.R")
# q()

if (FALSE) {

print("CS.mpp.R ..x.")
data.path <- " /home/schlather/rfm/data/"
dev <- 2;
cex<-3;cex.axis<-3;cex.lab<-3;cex.main<-3;cex.sub<-3;mar <- c(4,5,1,0.5)

lwd <- 5;

par(lwd=lwd)

######################################################################
#             Waelder & Stoyan, tree data set
######################################################################
name<-""; ## == dev=2 + readline()
SHOW<-TRUE
again <- FALSE;
nameEfct<-NULL;

WS <- as.matrix(read.table(paste(data.path,"waelder.stoyan.biondi.dat",sep=""),
                           header=T,skip=4))
plot(WS[,1],WS[,2])
nrow(WS) #149


######################################################################
#              Coulissenhieb
######################################################################
name<-"";
SHOW<-TRUE
again <- FALSE;
nameEfct<-NULL;


couliss <- read.table(paste(data.path,"Coulissenhieb.dat",sep=""), header=T,
                      skip=3)[,3:5]
Couliss <- as.matrix(couliss[!rowSums(is.na(couliss)),]) # remove NA's  
nrow(Couliss) #489
plot(Couliss[,c(1,2)])
Dev(T,dev=dev,"Couliss.diameter.ps",width=13,height=8)
par(lwd=5)
plot.with.circles(Couliss,factor=10,lwd=2,col=0,cex.axis=3)
Dev(F,dev=dev)

abline(20,0);abline(-60,0)
lines(c(-100,-100),c(-100,50))
lines(c(-20,-20),c(-100,50))
XCouliss<-Couliss[(Couliss[,1]> -100) & (Couliss[,1]< -20) &
                 (Couliss[,2]> -60) & (Couliss[,2]< 20), ]
points(XCouliss,pch="*")
nrow(XCouliss) #212

##plot(Couliss[,c(1,2)])
abline(20,0);
lines(c(-20,50),c(-20,-20))
lines(c(25,25),c(-20,20))
lines(c(-20,-20),c(-100,50))
YCouliss<-Couliss[(Couliss[,1]> -20) & (Couliss[,1]< 30) & (Couliss[,2]> -20) &
                 (Couliss[,2]< 20), ]
nrow(YCouliss) # 80

#plot(Couliss[,c(1,2)])
lines(c(-125,-85),c(35,35))
lines(c(-125,-125),c(-5,35))
lines(c(-125,-85),c(-5,-5))
lines(c(-85,-85),c(-5,35))
ZCouliss<-Couliss[(Couliss[,1]> -125) & (Couliss[,1]< -85) & (Couliss[,2]> -5) &
                 (Couliss[,2]<35), ]
nrow(ZCouliss) #56


######################################################################
#              Steigerwald, species not distinguished
######################################################################
name<-"";
SHOW<-FALSE
again <- FALSE;
nameEfct<-NULL;
showEfct <- TRUE

steiger <- read.table(paste(data.path,"Steigerwald.dat",sep=""),
                      header=T, skip=2)[, 3:6]
steiger <- steiger[!rowSums(is.na(steiger)),]
steiger[,3] <- steiger[,3]/2;
steiger <- steiger[steiger[,2]>min(steiger[,2]),] ## strange point
Steiger<- as.matrix(steiger[,c(1:3)])
nrow(Steiger) # 453
plot(Steiger[,c(1,2)], cex=cex, cex.axis=cex.axis, cex.lab=cex.lab,
     cex.main=cex.main, cex.sub=cex.sub)
Dev(T,dev=dev,"Steiger.diameter.ps",width=8,height=8,paper="special")
par(lwd=5)
plot.with.circles(Steiger,factor=10,lwd=2,col=0,cex.axis=3)
Dev(F,dev=dev)

# rotate the the location, then it is easier to get rid of the road
# visible in the data
phi <- +2/(2*pi)
TSteiger <- cbind( Steiger[,1]*cos(phi) - Steiger[,2] * sin(phi),
                 Steiger[,1]*sin(phi) + Steiger[,2] * cos(phi),Steiger[,3])
plot(TSteiger[,c(1,2)], cex=cex, cex.axis=cex.axis, cex.lab=cex.lab,
     cex.main=cex.main, cex.sub=cex.sub)

lines(c(-55,55),c(-15,-15))
lines(c(-55,55),c(45,45))
lines(c(-55,-55),c(-15,45))
lines(c(55,55),c(-15,45))
XTSteiger <- TSteiger[(TSteiger[,1]> -55) & (TSteiger[,1]<55 ) &
                     (TSteiger[,2]>-15 ) & (TSteiger[,2]<45 ), ]
nrow(XTSteiger) ## 229, point intensity in the nothern part: 229/(110*60)=0.03469

lines(c(-55,30),c(-20,-20))
lines(c(-55,30),c(-55,-55))
lines(c(-55,-55),c(-20,-55))
lines(c(30,30),c(-20,-55))
YTSteiger <- TSteiger[(TSteiger[,1]> -55) & (TSteiger[,1]< 30) &
                     (TSteiger[,2]> -55) & (TSteiger[,2]< -20), ]
nrow(YTSteiger) ## 142, point intensity in the southern part: 0.04773109

UpperSteiger<-TSteiger[TSteiger[,2]>-15,] #290
points(UpperSteiger[,c(1,2)],pch="*")
Dev(T,dev,ps="UpperSteiger.ps")
plot(UpperSteiger[,c(1,2)],pch="*", cex=cex, cex.axis=cex.axis, cex.lab=cex.lab,
     cex.main=cex.main, cex.sub=cex.sub)
xmin<--30; ymin<--10; deltaxy<-40
lines(c(xmin,xmin+deltaxy),c(ymin,ymin))
lines(c(xmin,xmin+deltaxy),c(ymin+deltaxy,ymin+deltaxy))
lines(c(xmin,xmin),c(ymin,ymin+deltaxy))
lines(c(xmin+deltaxy,xmin+deltaxy),c(ymin,ymin+deltaxy))
Dev(F,dev)
Dev(T,dev=dev,"UpperSteiger.diameter.ps",width=13,height=8)
plot.with.circles(UpperSteiger, factor=10, xlim=c(-60,70), ylim=c(-20,60),
                  lwd=2, col=0, cex.axis=3)
Dev(F,dev=dev)
ZTSteiger<-TSteiger[(TSteiger[,1]>xmin) & (TSteiger[,1]<xmin+deltaxy) &
                   (TSteiger[,2]>ymin) & (TSteiger[,2]<ymin+deltaxy), ]
print(nrow(ZTSteiger))

LowerSteiger<-TSteiger[TSteiger[,2]< -15,]
points(LowerSteiger[,c(1,2)],pch="*")



######################################################################
#              Steigerwald, species are distinguished
######################################################################

SteigerB <- Steiger[steiger[,4]=="B",]
SteigerE <- Steiger[steiger[,4]=="E",]
SteigerH <- Steiger[steiger[,4]=="H",]
print(c(nrow(Steiger),nrow(SteigerB),nrow(SteigerE),nrow(SteigerH)))

plot(Steiger[,1],Steiger[,2],col=0)
points(SteigerB[,1],SteigerB[,2],pch="B")
points(SteigerE[,1],SteigerE[,2],pch="E")
points(SteigerH[,1],SteigerH[,2],pch="H")

Usteiger <- steiger[  TSteiger[,2]>-15,4] 
USteigerB <- UpperSteiger[Usteiger=="B",]
USteigerE <- UpperSteiger[Usteiger=="E",]
USteigerH <- UpperSteiger[Usteiger=="H",]

Dev(T,dev=dev,"UpperSteiger.diameter.EBH.ps",width=13,height=8)
plot.with.circles(USteigerB, factor=10, xlim=c(-60,70), ylim=c(-20,60),
                  lwd=lwd, col=3)
plot.with.circles(USteigerE, factor=10, xlim=c(-60,70), ylim=c(-20,60),
                  lwd=lwd, col=1, new=TRUE)
Dev(F,dev=dev)

}
