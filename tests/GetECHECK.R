## source("GetECHECK.R")
if (EXTENDED.TESTING <- file.exists(f <- "CHECK.R")) source(f) else 
if (file.exists(f <- "~/R/MPP/MarkedPointProcess/tests/CHECK.R")) source(f) else
q()

#############



if (FALSE)
{
  load("delta.call")
  Delta(data=data,col=col,rep=rep,bin=bin,p=p,n=n,npoints=npoints,
        coordmodel=coordmodel)
  stop("")
}

cat("pid=", pid(), "hostname=", hostname(), "\n")


runif(1)
p <- 0.8
bin <- c(-1,seq(0,1.2,l=15))

npoints <- 20
z <- Delta(data=NULL,col=1,rep=1,bin=bin,p=p,n=1,npoints=npoints,
           coordmod="unif", cov.model="nugget")

npoints <- 20
z <- Delta(data=NULL,col=2,rep=2,bin=bin,p=p,n=1,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=1,rep=1,bin=bin,p=p,n=1,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=2,rep=1,bin=bin,p=p,n=1,npoints=npoints,coordmod="unif")

if (EXTENDED.TESTING) {
z <- Delta(data=NULL,col=1,rep=2,bin=bin,p=p,n=2,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=2,rep=2,bin=bin,p=p,n=1,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=2,rep=1,bin=bin,p=p,n=2,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=1,rep=2,bin=bin,p=p,n=2,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=2,rep=2,bin=bin,p=p,n=2,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=5,rep=1,bin=bin,p=p,n=1,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=1,rep=5,bin=bin,p=p,n=1,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=1,rep=1,bin=bin,p=p,n=5,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=5,rep=5,bin=bin,p=p,n=1,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=5,rep=1,bin=bin,p=p,n=5,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=5,rep=1,bin=bin,p=p,n=5,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=1,rep=5,bin=bin,p=p,n=5,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=5,rep=5,bin=bin,p=p,n=5,npoints=npoints,coordmod="unif")

z <- Delta(data=NULL,col=3,rep=5,bin=bin,p=p,n=7,npoints=npoints,coordmod="unif")
z <- Delta(data=NULL,col=6,rep=4,bin=bin,p=p,n=5,npoints=npoints,coordmod="unif")

lambda <- 20          
z <- Delta(data=NULL,col=1,rep=1,bin=bin,p=p,n=2,lambda=lambda,coordmod="Pois")
z <- Delta(data=NULL,col=1,rep=2,bin=bin,p=p,n=2,lambda=lambda,coordmod="Pois")
z <- Delta(data=NULL,col=2,rep=1,bin=bin,p=p,n=2,lambda=lambda,coordmod="Pois")
z <- Delta(data=NULL,col=1,rep=2,bin=bin,p=p,n=2,lambda=lambda,coordmod="Pois")
z <- Delta(data=NULL,col=2,rep=2,bin=bin,p=p,n=2,lambda=lambda,coordmod="Pois")

z <- Delta(data=NULL,col=1,rep=1,bin=bin,p=p,n=5,lambda=lambda,coordmod="Pois")
z <- Delta(data=NULL,col=5,rep=1,bin=bin,p=p,n=5,lambda=lambda,coordmod="Pois")

z <- Delta(data=NULL,col=1,rep=5,bin=bin,p=p,n=5,lambda=lambda,coordmod="Pois")
z <- Delta(data=NULL,col=5,rep=5,bin=bin,p=p,n=5,lambda=lambda,coordmod="Pois")
z <- Delta(data=NULL,col=3,rep=5,bin=bin,p=p,n=7,lambda=lambda,coordmod="Pois")
z <- Delta(data=NULL,col=6,rep=4,bin=bin,p=p,n=5,lambda=lambda,coordmod="Pois")

} # EXTENDED


if (FALSE)
{
  
  n <- 2;col <- 1;rep <- 2;
  p <- 0.8;bin <- c(-1,seq(0,1.2,l=15));npoints <- 20
  NEW  <-  FALSE
  if (NEW) {
    x  <- simulateMPP(npoints=npoints, coordrepet=n, window=c(0,1,0,1),
                       coordmod="unif", repetitions=rep, edgecorrection=0,
                       gauss.method="exponen", gauss.parameter=c(0,1,0,1))
    data  <-  list(); for (i in 1:length(x)) data  <-  c(data,x[[i]]) 
    save(x, data, bin, col, rep, p, file="GetECHECK.rda")
  } else load("GetECHECK.rda")
  z  <-  Delta(data=data,col=col,rep=rep,p=p,bin=bin)
  ##z  <-  mcfR(data=data,col=1,rep=2,p=0.8,bin=bin,analyseresult=NULL,getEresult=NULL)
  n <- 2;col <- 2;rep <- 1;
  z  <-  Delta(data=data,col=col,rep=rep,p=p,bin=bin,R=z$R) 
}
  

