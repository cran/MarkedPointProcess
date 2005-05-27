## source("MCTESTcheck.R")

if (EXTENDED.TESTING <- file.exists(f <- "CHECK.R")) source(f) else 
if (file.exists(f <- "~/R/MPP/MarkedPointProcess/tests/CHECK.R")) source(f) else
q()

#RFparameters(PrintLevel=2)
##EXTENDED.TESTING <- TRUE

bin <- c(-1,seq(0, 1.2, l=15));

cat("pid=", pid(), "hostname=", hostname(), "\n")

runif(1)

p <- c(0, 1, 0, 1)

pts <- 20
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG,Barnar=TRUE)
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG,Barnar=TRUE)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG,Barnar=TRUE)

z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)

if (EXTENDED.TESTING) {
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)


lambda <- 20
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,lam=lambda,coordm="Pois",par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,lam=lambda,coordm="Pois",par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,lam=lambda,coordm="Pois",par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,lam=lambda,coordm="Pois",par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,lam=lambda,coordm="Pois",par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,lam=lambda,coordm="Pois",par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,lam=lambda,coordm="Pois",par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,coordmodel=2,par=p)

pts <- 50
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 15
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 10
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 5
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 50
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 15
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 10
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 5
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 50
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 15
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 10
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 5
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 50
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 15
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 10
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 5
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 50
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 15
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 10
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 5
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 50
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 15
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 10
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)

pts <- 5
z <- DeltaMC(rep=1,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=1,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=1,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=2,bin=bin,n=2,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=7,bin=bin,n=5,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
z <- DeltaMC(rep=5,bin=bin,n=7,MCrep=99,npoints=pts,par=p,DEB=DEBUG)
#DeltaMC(rep=30,bin=bin,n=30,MCrep=99,npoints=pts,par=p)


} # EXTENDED





