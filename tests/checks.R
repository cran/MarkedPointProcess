## source("checks.R")

# Here, some checks are performed that are not covered by GetECHECK.R and MCTESTcheck.R,
# i.e., checks if the simulations of the random fields are OK
# (In GetECHECK.R and MCTESTcheck.R, the calculation of the Efunction and the
# statistics are checked by comparing the results obtained by the C programme and
# the (very slow) R implementation )


##################################################################################
## Rechecked part

##################################################################################
## is covariance structure given correctly?
## uncomment #define debug in auxiliary first

if (EXTENDED.TESTING <- file.exists(f <- "CHECK.R")) source(f)


runif(1)

RFparameters(Print=1)

individ <- 50
radius <- 0.1
xlim <- c(0,1)
ylim <- c(0,1)
bin <- c(-1,seq(0,0.7,l=15))
sill <- 0.8


simu <- function(repet, MCrep) {
  if (repet==1) stop(" at least 2")
  x <- simulateMPP(npoints=individ ,coordmodel="unif", window=c(xlim,ylim),
                    edgecorrection = radius, repetitions=1, coordrepet=repet,
                    model=
                   list("+",
                      list("$", var=sill, scale=0.1, list("spherical")),
                      list("random coin",
                           p=c(1, radius, -sqrt((1-sill) *
                             diff(xlim)*diff(ylim)/(individ*pi*radius^2))) ))
                   )
  res <- list();
  r <- matrix(nrow=length(MCrep),ncol=5)
  RFparameters(Print=1)
  for (i in 1:length(MCrep)) {
    cat(i, "")
    r[i,] <- unix.time(res[[i]] <-
                       rfm.test(data=x, MCrep=MCrep[i], bin=bin, pvalue=NULL))
  }
  return(list(r=r, res=res, x=x))
}
#
if (FALSE) {
#options(warn=2)
res0 <- simu(rep=10, MC=c(99, 199, 499))
stop("")
}

res0 <- simu(rep=3, MC=c(99, 199))

if (EXTENDED.TESTING) {
res1 <- simu(rep=100, MC=c(99, 199, 499))
res2 <- simu(rep=500, MC=c(99, 199, 499))
res3 <- simu(rep=1000, MC=c(99, 199, 499))
res4 <- simu(rep=2000, MC=c(99, 199, 499))
res5 <- simu(rep=5000, MC=c(99, 199, 499))

res11 <- simu(rep=100, MC=c(99, 199, 499))
res12 <- simu(rep=500, MC=c(99, 199, 499))
res13 <- simu(rep=1000, MC=c(99, 199, 499))
res14 <- simu(rep=2000, MC=c(99, 199, 499))
res15 <- simu(rep=5000, MC=c(99, 199, 499))


res21 <- simu(rep=100, MC=c(99, 199, 499))
res22 <- simu(rep=500, MC=c(99, 199, 499))
res23 <- simu(rep=1000, MC=c(99, 199, 499))
res24 <- simu(rep=2000, MC=c(99, 199, 499))
res25 <- simu(rep=5000, MC=c(99, 199, 499))

k <- 1
x <- seq(0,0,l=100); n <- 1;
for (i in 1:100) {
  if  (res1[[2]][[1]]$E[i,k]>0) {
    for (j in 1: res1[[2]][[1]]$E[i,k]) {
      x[n] <- i;
      n <- n+1
    }
  }
}
mean(x)
sqrt(var(x))

nn  <- c(100,500,1000,2000,5000)
for (k in 1:1) {
  z <- NULL
  for (i in 1:3) {
    result <-  eval(parse(text=paste("res",k,"[[2]][[i]]",sep="")))
    m <-  colSums(as.double(1:100) * result$E)
    m2 <- colSums(as.double(1:100)^2 * result$E)
    print(m)
    print(m2)
    s   <- sqrt((m2-m * m/nn[k])/(nn[k]-1)/nn[k])
    m <- m / nn[k]
    z <- cbind(z,m,s=nice(s,dec=3),"    ")
  }
  print(z,quote=FALSE)
  readline()
}

k <- 1;cbind(res4[[2]][[1]]$E[,k],res4[[2]][[2]]$E[,k],res4[[2]][[3]]$E[,k])
k <- 1;cbind(res3[[2]][[1]]$E[,k],res4[[2]][[1]]$E[,k]/2)
plot(cumsum(res3[[2]][[1]]$E[,k]),cumsum(res4[[2]][[1]]$E[,k]/2))


E <- matrix(nrow=20,ncol=NUMBERTESTS)
for (i in 1:20)
 E[i,] <-  rfm.test(x$coord,x$data, MC=99, bin=bin)$E

t(E)

E2 <- matrix(nrow=20,ncol=NUMBERTESTS)
for (i in 1:20)
 E2[i,] <-  rfm.test(x$coord,x$data, MC=199, bin=bin)$E

sqrt(apply(E,2,var))
sqrt(apply(E2/2,2,var))

colMeans(E)
colsMenas(E2/2)


} # EXTENDED







######################################################################
######################################################################
######################################################################
