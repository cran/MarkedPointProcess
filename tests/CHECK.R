
## source("CHECK.R")

## only works if Epos=1,2,3,4,...,100 (standard!)
## and MCrep=99

epsilonvariance <- 1E-7 ##check with MPPanalyse.cc -- same values of the variables!!

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

DEBUG <- FALSE; #DEBUG <- TRUE
AGAIN <- DEBUG



runif(1)
dummy <- .C("GetmppParameters", lnorms=integer(1), weights=integer(1),
            tests=integer(1), mppmaxchar=integer(1), modelnr=integer(1),
            PACKAGE="MarkedPointProcess", DUP=FALSE)
.mpp.l.norms <- dummy$lnorms
.mpp.weights <- dummy$weights
.mpp.tests <- dummy$tests
.mpp.maxchar <- dummy$mppmaxchar
.mpp.models <- dummy$modelnr


teststatistic <- function(E, Ebin, sigma, sd, col, n, rep, p) {
  weight <- matrix(nrow=nrow(Ebin), ncol=.mpp.weights)
  result <- matrix(NA, nrow=.mpp.tests, ncol=ncol(E))
  Ebin <- as.matrix(Ebin)
  E <- as.matrix(E)
  repeatindex <- ((0:(rep-1)) * (col*n*n))
  sdrepindex  <- ((0:(rep-1)) * (col*n))
  zaehler <- 0 ##
  for (i in 1:n) {
    for (j in 1:n) {
      for (cc in 1:col) {
        zaehler <-  zaehler + 1 ##
        majorindex <- (i-1)*n*col + (j-1)*col + cc;
        zeroindex <- (i-1)*(n+1)*col + cc ## E[0]

        SDa <- sd[(i-1)*col + cc + sdrepindex]        
        EE <- t(E[, majorindex + repeatindex, drop=FALSE])## achtung, bins quer--
        EEbin <- Ebin[, majorindex, drop=FALSE]
        weight[,1] <- 1
        weight[,5] <- EEbin; 
        weight[,6] <- sqrt(weight[,5])
        weight[,4] <-
          c(0, 1/cumsum(c(sqrt(Ebin[1,zeroindex]), weight[-1,6])[-nrow(Ebin)]))
        weight[,2] <-
          c(0, 1/cumsum(c(Ebin[1,zeroindex], weight[-1,5])[-nrow(Ebin)]))
        weight[,3] <- sqrt(weight[,2])
        weight[,7] <- 1 ## dummy value so that apply doesn't do wrong,
        ##                 true value is set later on!
        weight <- t(t(weight)/apply(weight[-1,], 2, sum))
        ## if ((testname=="Stest") && (zaehler==-20))  print(weight)
        EE0 <-  abs(EE-E[1,repeatindex+zeroindex])  ## achtung, bins quer ---
        tEE0 <- t(EE0)                              ## achtung, bins laengs |
        ## weight == 7
        sumsigma <- apply(1/sigma[-1,majorindex + repeatindex,drop=FALSE],2,
                          sum, na.rm=TRUE)
        
        tSIGMA0 <- t(t(1/sigma[, majorindex+repeatindex, drop=FALSE]) / sumsigma)
               
        # "ROBUST":
        quantile<-function(x,probs,na.rm=NULL){
          x <- x[is.finite(x)]
          (x[order(x)])[probs*length(x)+1];
        }
        q <- apply(cbind(0, EE0[, -1, drop=FALSE]), 1, quantile, na.rm=TRUE,
                   probs=p)

        if ((testname=="Etest") && (zaehler==-33))
          {cat("\n",zaehler,"*** q.cc: ",q,"\n")}
        
        if (DEBUG) print(c("QQQ",q))
        if (DEBUG) print(EE0)
        ind <- !is.na(EE0) & EE0 < q
        
        rq <- (2 * q) * EE0 - q^2
        rq[ind] <- EE0[ind]^2
        assign("ll", list(ind, EE0, rq, q), envir=.GlobalEnv)
        resi <- 1
 
        for (w in 1:.mpp.weights) {
          if (w<7) { gewicht <- weight[,w] }
          else {
            gewicht <- tSIGMA0
          }
          result[resi, majorindex + repeatindex] <-
            apply(gewicht * tEE0, 2, max, na.rm=T)
          resi <- resi+1
          if (DEBUG) { print("w",w); print(gewicht*tEE0) }
          result[resi, majorindex + repeatindex] <-
            apply(gewicht * tEE0^2, 2, sum, na.rm=T)
          resi <- resi+1
          result[resi, majorindex + repeatindex] <-
            apply(gewicht * tEE0, 2, sum, na.rm=T)
          resi <- resi+1
          result[resi, majorindex + repeatindex] <-
            apply(gewicht * t(rq), 2, sum, na.rm=T)
##if (w==7) print(result[resi, majorindex + repeatindex])
          resi <- resi+1
          for (s in (1:1)/10) {
             r <- EE0
             ax <- SDa*s
             bx <- 0.25 / ax
             ind <- !is.na(r) & r >= ax
             r[ind] <- (bx * (r + ax)^2) [ind]
             result[resi, majorindex + repeatindex] <-
               apply(gewicht * t(r), 2, sum, na.rm=T)
             resi <- resi+1;           
          }
        }        
      } #c
    } #j
  } #i
  result[resi,] <- apply(as.matrix(apply(E, 2, range, na.rm=T)), 2, diff)
  return(as.matrix(result))
}


mcfR <- function(data, col, rep, p, bin, analyseresult=NULL, getEresult=NULL) {
  for (i in 1:length(data)) data[[i]] <- as.matrix(data[[i]])
  n <- as.integer(length(data)/2) ## length of list [coord,data] == #species
  nb <- length(bin)-1
  coltri <- col * (col+1) /2
  estSd <- estE <- double(rep * col * n) ## so zero by default
  estSDsd <- estSDVar <- double(rep * coltri * n) 
  index <- as.vector(t(t(matrix(1:col,nrow=col,ncol=rep)) + (0:(rep-1))*n*col))
  ## 1, 1+n, (1+n)+(n-1), ..., ()+2

  var.index <- as.vector(t(t(matrix( if (col==1) 1 else c(1, 1+cumsum(col:2)),
                                    nrow=col, ncol=rep)) +
                           (0:(rep-1))*n*coltri))
  for (i in 1:n) {
    estE[ index + (i-1)*col ] <- apply(data[[2*i]], 2, mean)
    estSd[ index + (i-1)*col ] <- sqrt(apply(data[[2*i]], 2, var))
  }
  for (i in 1:n) {
    classes <- as.integer(sqrt(nrow(data[[2*i]])))    
    classsize <- as.integer(nrow(data[[2*i]]) / classes)
    classbounds <- (1:classes)*classsize
    deltaclass <- nrow(data[[2*i]]) - classbounds[length(classbounds)]
    if (deltaclass > 0) classbounds <- classbounds +
      c(cumsum(seq(1, 1, l=deltaclass)),
        seq(deltaclass, deltaclass, l = classes - deltaclass))
    classbounds <- c(0, classbounds)
    x <- NULL
    for (j in 1:classes) {
      x <- rbind(x, apply(data[[2*i]][(classbounds[j]+1):classbounds[j+1],,
                                      drop=FALSE],2,var) )
    }
    ## estimated standard deviation of the standard deviation
    estSDsd[var.index + (i-1)*coltri] <- sqrt(apply(sqrt(x),2,var)/classes)
    ## estimated standard deviation of the variance
    estSDVar[var.index + (i-1)*coltri] <- sqrt(apply(x,2,var)/classes)
  }
  
  E   <- ESIGMA <- matrix(nrow=nb,ncol=n*n*rep*col)
  VAR <-VARSIGMA <- V22 <- V21<- V12<- V11<-
    matrix(nrow=nb, ncol = n * n * rep * col * (col+1)/2)
  KMM <- matrix(nrow=nb, ncol = rep * (n * col) * (n * col + 1) / 2)
  GAM <- matrix(nrow=nb, ncol = rep * (n * col) * (n * col + 1) / 2)
  EBIN<- matrix(nrow=nb, ncol = n * n * col)
  VARBIN <- matrix(nrow=nb, ncol = n * n * col * (col + 1) / 2)
  KMMBIN <- matrix(nrow=nb, ncol = (n * col) * (n * col + 1) / 2)
  GAMBIN <- matrix(nrow=nb, ncol = (n * col) * (n * col + 1) / 2)

  repeatindex <- (0:(rep-1)) * col;
  erepindex   <- (0:(rep-1)) * col * n * n;
  varrepindex <- (0:(rep-1)) * coltri * n * n;
  result <- matrix(0, nrow=2, ncol = n * n * col * rep)  
  for (s1 in 1:n) {
    n1 <- nrow(data[[2*s1-1]])
    for (s2 in 1:n) {
      n2 <- nrow(data[[2*s2-1]])
      x1 <- matrix(1:n1,nrow=n1,ncol=n2)
      distances <- sqrt(outer(data[[2*s1-1]][,1], data[[2*s2-1]][,1], "-")^2 +
                        outer(data[[2*s1-1]][,2], data[[2*s2-1]][,2], "-")^2)
      for (b in 1:nb) {
        index <- x1[(distances>bin[b]) & (distances<=bin[b+1])]
        for (c1 in 1:col) {
          if (length(index>0)) {
            E[b, c1 + (s2-1)*col + (s1-1)*n*col + erepindex] <-
              apply(data[[2*s1]][index,c1+repeatindex,drop=FALSE],2,mean)
            ESIGMA[b, c1 + (s2-1)*col + (s1-1)*n*col + erepindex] <- 
              sqrt(apply(data[[2*s1]][index,c1+repeatindex,drop=FALSE],2,var))
            ## for replacement of NA by Inf see below
          }
          EBIN[b, c1 + (s2-1)*col + (s1-1)*n*col] <- length(index)
          ##                                                     NA not allowed !
        }
        k <- 1
        for (c1 in 1:col) {
          for (c2 in c1:col) {
            if (length(index>1)) {
              for (r in 1:rep) {
                #(c1=1,c2=1), (c1=1,c2=2), (c1=2, c2=2)
                VAR[b, k + (s2-1)*coltri + (s1-1)*n*coltri + varrepindex[r]] <-
                  cov(data[[2*s1]][index,c1 + repeatindex[r]],
                      data[[2*s1]][index,c2 + repeatindex[r]])
              }
              if (length(index)>2) { for (r in 1:rep) {               
                VARSIGMA[b, k +(s2-1)*coltri +(s1-1)*n*coltri + varrepindex[r]]<-
                  sqrt(sum(((data[[2*s1]][index, c1 + repeatindex[r]] -
                             mean(data[[2*s1]][index,c1+repeatindex[r]])) 
                            * (data[[2*s1]][index, c2 + repeatindex[r]] -
                               mean(data[[2*s1]][index, c2 + repeatindex[r]]))
                            - VAR[b, k+(s2-1)*coltri+(s1-1)*n*coltri
                                  + varrepindex[r]]
                            )^2 ) / (length(index)-2))
              }} # length(index)>2, for r
            }
            VARBIN[b, k + (s2-1)*coltri + (s1-1)*n*coltri] <- length(index)
            k <- k+1;
          }
        }
      } # b
    } # s2
  } #s1
  SQ <- sqrt(abs(VAR)) * sign(VAR)
  SQBIN <- VARBIN
  
  k<-1;
  kmmrepindex <- (0:(rep - 1)) * (n * col) * (n * col + 1) /2;
  for (i in (0:(n * col - 1))) {
    s1 <- as.integer(i / col);
    c1 <- (i - s1 * col) + 1
    s1 <- s1 + 1
    n1 <- nrow(data[[2 * s1 - 1]])
    for (j in (i:(n * col-1))) {
      s2 <- as.integer(j / col);
      c2 <- (j - s2 * col) + 1
      s2 <- s2 +1
      n2 <- nrow(data[[2 * s2 - 1]])
      x1 <- matrix(1:n1, nrow=n1, ncol=n2)
      x2 <- t(matrix(1:n2, nrow=n2, ncol=n1))
      distances <- sqrt(outer(data[[2*s1-1]][,1], data[[2*s2-1]][,1], "-")^2 +
                        outer(data[[2*s1-1]][,2], data[[2*s2-1]][,2], "-")^2)
      for (b in 1:nb) {
        index1 <- x1[(distances > bin[b]) & (distances <= bin[b + 1])]
        index2 <- x2[(distances > bin[b]) & (distances <= bin[b + 1])]
        GAMBIN[b, k] <- KMMBIN[b, k] <- length(index1)
        if (length(index1 > 0)) {
          KMM[b, k + kmmrepindex] <-
            apply(data[[2 * s1]][index1, c1 + repeatindex, drop=FALSE] *
                  data[[2 * s2]][index2, c2 + repeatindex, drop=FALSE], 2, mean)
          GAM[b, k + kmmrepindex] <-
            apply((data[[2 * s1]][index1, c1 + repeatindex, drop=FALSE] -
                   data[[2 * s2]][index2, c1 + repeatindex, drop=FALSE]) *
                  (data[[2 * s1]][index1, c2 + repeatindex, drop=FALSE] -
                   data[[2 * s2]][index2, c2 + repeatindex, drop=FALSE]),
                  2, mean) / 2
        }
      }
      k <- k + 1
    }
  }
  et <- vt <- st <- et.rank <- vt.rank <- st.rank <- NULL

  assign("testname", "Etest", envir=.GlobalEnv)

  ## in case the variables (take discrete values and thus can) have all the
  ## same
  ESIGMA[ESIGMA<epsilonvariance] <- epsilonvariance
  
  et <- teststatistic(E, EBIN, ESIGMA, sd=estSd, col=col, n=n, rep=rep, p=p)
  if (ncol(et) > 2) et.rank <- apply(et[,1] >= et, 1, sum)

  assign("testname", "Vtest", envir=.GlobalEnv)
  VARSIGMA[VARSIGMA<epsilonvariance] <- epsilonvariance
  vt <- teststatistic(VAR, VARBIN, VARSIGMA, sd=estSDVar, col=coltri, n=n,
                      rep=rep, p=p)
  vt.rank <- matrix(vt, nrow=nrow(vt) * coltri * n * n)
  if (ncol(vt) > 2) vt.rank <- apply(vt[,1] >= vt, 1, sum)
  
  assign("testname", "Stest", envir=.GlobalEnv)
  st <- teststatistic(SQ, SQBIN, sqrt(VARSIGMA), sd=estSDsd, col=coltri, n=n,
                      rep=rep, p=p)
  if (ncol(st) > 2) st.rank <- apply(st[, 1] >= st, 1, sum)
  l <- list(E=E, Ebin=EBIN, ETest=et, ET=et.rank,
            VAR=VAR, VARbin=VARBIN, VARTest=vt, VT=vt.rank,
            SQ=SQ, SQbin=SQBIN, SQTest=st, SQT=st.rank,
            KMM=KMM, KMMbin=KMMBIN, GAM=GAM, GAMbin=GAMBIN)
  return(l)
}


Delta <- function(data, col, rep, bin, p=0.8, n=NULL, npoints=NULL,
                  lambda=npoints/2, edgecorrection=0,
                  coordmodel="given", window=c(0, 1, 0, 1),
                  R=NULL, C=NULL, cov.model="exponential"){
  ## test
  ##
  ## the function and statistics for data are calculated and the results of the
  ## C and the (slow) R implementation are compared
  ##
  ## data : null : data are first simulated
  ##        list(coord,data,coord,data,...) expected, otherwise
  ## col : col-variate data
  ## rep : how often repeated and analysed simulataneously
  ## bin : "variogram" bin
  ## ADDITIONAL PARAMETERS FOR SIMULATION
  ## n : number of species
  ## npoints : (maximal) number of locations, see simulateMPP
  ## coordmodel, see simulateMPP
  ## R : not NULL : Delta$R expected
  ## C : not NULL : Delta$C expected
  ##
  ## returns R,C,data
    
  data; npoints; lambda; R; C;
  cat("col=",col,"rep=",rep,"n=",n,"npoints=",npoints,"coordmodel=",coordmodel,
      "\n")
  
  save(file="delta.call",
       .Random.seed, data, col, rep, bin, p, n, npoints,
       coordmodel, window, lambda, edgecorrection, R, C)

  pos.variance <- TRUE
  if (is.null(data)) {
    R <- C <- NULL;
    x <- simulateMPP(npoints=npoints, coordrepet=n, window=window,
                      lambda=lambda, edgecorrection=edgecorrection,
                      coordmodel=coordmodel, repetitions=rep*col,
                      model = list(list(model=cov.model, var=1, scale=1)) )
  norm <- function(x) {qnorm((rank(x)-0.5)/length(x))}
   if (n==1) {
      if (any(apply(x$data,2,var))==0) pos.variance <- FALSE
      else  x$data <- apply(as.matrix(x$data),2,norm)
      data <- x
     }
    else {
      data <- list();
      for (i in 1:length(x)) {
        if (any(apply(x[[i]]$data,2,var))==0) pos.variance <- FALSE
        else x[[i]]$data <- apply(as.matrix(x[[i]]$data),2,norm)
        data <- c(data, x[[i]])
      }
    }
  }

   if (!pos.variance) {
    print(data)
    print("^^^^^^^^^^^^^^^^^^^^^   !pos variance")
    stop("")
    return();
  }
           
  ##  print("got Data")
  if (is.null(C))
    C <- mpp.characteristics(bin=bin, rep=rep, p=p, summarize=FALSE,
                             normalize=FALSE, data=data, staticchoice=TRUE)
  
  if (is.null(R)) R <- mcfR(data=data, col=col, rep=rep, p=p, bin=bin)

  for (name in names(C)) {
    eval(parse(text=paste("RX <- R$",name)))
    eval(parse(text=paste("CX <- C$",name)))
    ##if (name=="VAR") stop("")
    if (!is.null(RX)) {
      if (dim(RX)[1]==38) {
        idx <- c(5, 10, 15, 20, 25, 30, 35)
        id <- c(1:4, 6:9, 11:14, 16:19, 21:24, 26:29, 31:34)
        RX[37:38,] <- 0;
        CX[37:38,] <- 0
      } else {
        idx <- NULL
        id <-  1:nrow(RX)
      }
      z <- RX - CX
      if (name=="SQ") thr <- 1e-7  else thr <- 1e-13
      if (name=="SQTest") {
        rel.thrX <- 5e-3
        rel.thr <- 5e-3
      } else {
        rel.thrX <- 1e-7
        rel.thr <- 1e-7
      }
      s <- abs(z/RX) * ((abs(RX)>thr) | (abs(CX)>thr))
      nas <- sum(abs(is.finite(RX) - is.finite(CX)))
      if   ((any(s[id, ] > rel.thr, na.rm=TRUE)) ||
            (any(s[idx, ] > rel.thrX, na.rm=TRUE)) ||
            (nas>0)) {
        cat("\n$$$$$$$$ col=", col, "rep=", rep, "bin=", bin, "p=", p, "n=", n,
            "npoints=", npoints ,"coordmodel=", coordmodel,"\n")
        dd <- (
               1 * (is.finite(RX+CX) & ((abs(RX)>thr) | (abs(CX)>thr)) &
                    (abs(z/RX)>rel.thr))
               + 2 * (is.finite(RX) & !is.finite(CX))
               + 4 * (!is.finite(RX) & is.finite(CX))
               + 8 * (!is.finite(RX+CX) & is.na(RX) & !is.na(CX))
               +16 * (!is.finite(RX+CX) & !is.na(RX) & is.na(CX))
               +32 * (!is.finite(RX+CX) & !is.na(RX) & !is.na(CX)
                      & (RX==Inf) & CX==-Inf)
               +64 * (!is.finite(RX+CX) & !is.na(RX) & !is.na(CX)
                      & (RX==-Inf) & CX==Inf)
               )
        slct <- dd>0
        print(cbind(matrix(1:nrow(dd),nrow=nrow(dd), ncol=ncol(dd))[slct],
                    t(matrix(1:ncol(dd),nrow=ncol(dd), ncol=nrow(dd)))[slct],
                    dd[slct],
                    RX[slct],
                    CX[slct],
                    z[slct],
                    s[slct]
                    ))
        print(nas)
      }
      if (FALSE) {
        warn1 <- 2e-8
        warn2 <- 1e-14
        if (any(s>warn1, na.rm=TRUE) || (nas>0))
          cat("*********** WARNING  *******\n")
        else if (any(s>warn2, na.rm=TRUE))  cat("+++++++++++ WARNING  +++++++\n")
        cat(name,
            ": m1=",mean(s[s>warn1],na.rm=TRUE),"[",sum(s>warn1,na.rm=TRUE),"]",
            " m2=",mean(s[s>warn2],na.rm=TRUE),"[",sum(s>warn2,na.rm=TRUE),"]",
            " max=", max(s,na.rm=TRUE),
            " NA/notNA=", nas, "\n")
      }
    } else {
     # cat(name,"is unknown!\n")
    }
  }
  return(list(R=R, C=C, data=data))
}

DeltaMC <- function(data=NULL, rep, bin, n, MCrep,
                    npoints=NULL, coordmodel="unif",
                    window=c(0, 1, 0, 1), edgecorrection=0, lambda=npoints/2,
                    DEBUG=FALSE, again=AGAIN,
                    model="exponential",
                    param=c(mean=0, variance=NA, nugget=0, scale=NA),
                    sill=NA, Barnard=FALSE){
  
  ## NOTE: IN ORDER TO USE THIS FUNCTION COMPILE mpp.cc WITH #define DEBUG 1
  ## ********************************************************************
  ## test
  ## the function and statistics for data are calculated and the results of the
  ## C and the (slow) R implementation are compared
  ##
  ## data : null : data are first simulated
  ##        list(coord,data,coord,data,...) expected, otherwise
  ## rep : how often repeated and analysed simulataneously
  ## bin : "variogram" bin
  ## n : n-variate data
  ## ADDITIONAL PARAMETERS FOR SIMULATION
  ## npoints : (maximal) number of locations, see simulateMPP
  ## coordmodel, see simulateMPP
  ## R : not NULL : Delta$R expected
  ## C : not NULL : Delta$C expected
  ##
  ## returns R,C,data

  RFparameters(Storing=TRUE)
  seedfile <- paste("DeltaMC", rep, n, npoints, lambda, "rda",sep=".")
  if (again) {
    load(seedfile)
    assign(".Random.seed", seed, envir=.GlobalEnv)
  } else {
    seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    save(file=seedfile, seed)
  }
    
  p <- 0.7 ## values as defined in MCtest

  print(param)
  
  if (is.vector(param))  param <- matrix(param, ncol=n, nrow=length(param));
#  if (ncol(param) < rep) param <- matrix(param, nrow=nrow(param), ncol=rep)

  if (DEBUG) cat("simulate...")
  if (is.null(data)) {
    data <-
      simulateMPP(npoints=npoints, coordrepet=n, window=window,
                   edgecorrection=edgecorrection, lambda=lambda,
                   coordmodel=coordmodel, repetitions=rep,
                   model = list(list(model=model,  var=1, scale=1)) )
    if (n==1) data <- list(data)
  }
  ## rfm.test(data=data, normalize=FALSE, MCrep=MCrep,
  ##         MCmodel=model, MCparam=param, sill=sill, bin=biin,
  ##         use.naturalscaling=TRUE)
  
  lEbinM1 <- MCrep + 1
  Etest <- integer(lEbinM1 * .mpp.tests)
  VARtest <- integer(lEbinM1 * .mpp.tests)
  SQtest <- integer(lEbinM1 * .mpp.tests) ###
  MAXtest <- integer(lEbinM1 * .mpp.nr.maxtests) ###

  ## arrays to store places with differences
  ET <- matrix(0, nrow=lEbinM1, ncol=.mpp.tests)
  VT <- matrix(0, nrow=lEbinM1, ncol=.mpp.tests)
  SQT <- matrix(0, nrow=lEbinM1, ncol=.mpp.tests)
  ##MT = matrix(MAXtest,nrow=length(Ebin)-1)
  additive <- as.integer(1)
  error <- integer(1)

  EPSILON <- 0.000000001
  intestbin <-  1 : (MCrep + 1)

  m <- list()
  if (DEBUG) cat("fitvario...")
  for (r in 1:length(data)) {
    Data <- as.matrix(data[[r]]$data)
    for(i in 1:ncol(Data)) {
      if (DEBUG) cat(i, "")
      
      simuresult <- if (Barnard) {
        Data[sample(nrow(Data), size=nrow(Data) * MCrep, replace=TRUE), i]
      } else {
        GaussRF(x=data[[r]]$coord, grid=FALSE,
                model=fitvario(x=data[[r]]$coord, data=Data[,i], model=model,
                  param=param[,r], sill=sill)$variogram$ml,
                register=2, n=MCrep)
      }
      storage.mode(simuresult) <- "double"
      
      if (DEBUG) cat("start mctest...")

      .C("MCtest",
         as.integer(MCrep),	
         as.double(data[[r]]$coord),
         as.double(Data[,i]),
         as.integer(nrow(Data)),
         as.integer(ncol(data[[r]]$coord)),
         simuresult, 
         as.integer(1 + 4 * DEBUG), # PrintLevel
         as.double(bin),
         as.integer(length(bin)-1),
#         as.double(Ebin),
#         as.integer(length(Ebin)-1),
         Etest,
         VARtest,
         SQtest,  ###
         as.integer(.mpp.maxtests),
         as.integer(.mpp.nr.maxtests),
         MAXtest,
         error,
         additive,
         as.integer(1), # copy simulation result if barnard & static choice in mcf_iinternal
         PACKAGE="MarkedPointProcess", DUP=FALSE
         )
      if (error) stop(paste("error #",error,"has occured"))
      stopifnot(!is.null(simuresult))
      
      Etest <- matrix(Etest, ncol=.mpp.tests)
      VARtest <- matrix(VARtest, ncol=.mpp.tests)
      SQtest <- matrix(SQtest, ncol=.mpp.tests)
      simuresult <- matrix(simuresult, nrow=nrow(Data))    
      
      if (DEBUG) cat("testing...")
      if (DEBUG) cat(i, "")
      daten <- cbind(Data[,i], simuresult)
      x <- mcfR(data=list(data[[r]]$coord,
                  daten), #+ rnorm(length(daten), 0, 0.0001)), # force errors
                  col=1, rep=MCrep+1, p=p, bin=bin)
      for (j in 1:.mpp.tests) {
        ET[x$ET[j], j]   <- ET[x$ET[j], j] + 1
        VT[x$VT[j], j]   <- VT[x$VT[j], j] + 1
        SQT[x$SQT[j], j] <- SQT[x$SQT[j], j] + 1
      }
    }
  }

#  print("FINAL")
  if (DEBUG) {
    print(ET-Etest)
    print(VT-VARtest)
    print(SQT-SQtest)
    print("simuresult")
    print(simuresult)
  }

#  print(ET) #12
#  print(Etest) #13

  print(match.call())
  print("summaries")
  sE <- sum(abs(ET-Etest)[, 1:(.mpp.tests-2)], na.rm=TRUE)
  sV <- sum(abs(VT-VARtest)[, 1:(.mpp.tests-2)], na.rm=TRUE)
  sS <- sum(abs(SQT-SQtest)[, 1:(.mpp.tests-2)], na.rm=TRUE)
  print(sE)
  print(sV)
  print(sS)
  if (sE + sV + sS>0) {
      l <- list(data=data, simuresult=simuresult, rep=rep, bin=bin, n=n,
                MCrep=MCrep, npoints=npoints, coordmodel=coordmodel, m=m,
                window=window, edgecorrection=edgecorrection,
                lambda=lambda, ET=ET, Etest=Etest,
                VT=VT, VARtest=VARtest, SQT=SQT, SQtest=SQtest,
                E=x$E, VAR=x$VAR, SQ=x$SQ, seed=seed
               )
      save(l,x,file="error.x") #   load("error.x"); sum(abs(l$VT-l$VARtest)[,1:(.mpp.tests-2)],na.rm=TRUE); l$VT-l$VARtest
      ERROR    
      # if (EXTENDED.TESTING) ERROR else warning("ERROR: Differences observed!")
   }
}
