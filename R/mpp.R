

############# WORK TO DO ##################################################
## check why V fucntion goes crazy -- as mentioned in the SRD-Paper
###########################################################################

#
.ENV <- environment(NULL)
#.ENV <- baseenv()
#.ENV <- .GlobalEnv

## basic initializations in R
assign(".mpp.maxtests",
       as.integer(t(as.matrix(expand.grid(c(12,20,14,18,23,27),
                                          c(12,13,21,22,30,31))))), env=.ENV)
assign(".mpp.nr.maxtests", as.integer(length(.mpp.maxtests) / 2), env=.ENV)
assign(".mpp.digits", 3, env=.ENV)
assign(".mpp.lpnames", c("max", "l2", "l1", "robust", "anti"), env=.ENV)
assign(".mpp.weightnames", c("const", "1/sum#", "sqrt(1/sum#)", "1/sumsqrt#",
                             "#", "sqrt#", "1/sd") , env=.ENV)
assign(".mpp.extranames", c("range", "no.bin.sq", "no.bin.abs"), env=.ENV)

.onLoad <- function (lib, pkg) {
  if (file.exists("/home/schlather/bef/x")) {
    ## to-do list -- since my surname is rare, the message should 
    ## appear only on computers I have a login
    cat("To-Do List\n==========\n")
  }
 
  ##  dummy <- .C("GetmppParameters", lnorms=integer(1), weights=integer(1),
  ##               tests=integer(1), mppmaxchar=integer(1), modelnr=integer(1),
  ##               PACKAGE="MarkedPointProcess", DUP=FALSE)
  ## .mpp.l.norms <- dummy$lnorms
  ## .mpp.weights <- dummy$weights
  ## .mpp.tests <- dummy$tests
  ## .mpp.maxchar <- dummy$mppmaxchar
  ## .mpp.models <- dummy$modelnr
  
  ## see also tests/CHECK.R for GetmppParameters !!
  ##
  #  assign("..ENV", .ENV, envir=.ENV)
}

get.mpp.names <- function() {
  dummy <- .C("GetmppParameters", lnorms=integer(1), weights=integer(1),
                     tests=integer(1), mppmaxchar=integer(1), modelnr=integer(1),
                     PACKAGE="MarkedPointProcess", DUP=FALSE)
  .mpp.models <- dummy$modelnr
  .mpp.maxchar <- dummy$mppmaxchar
  l <- NULL;
  for (i in 1:.mpp.models) {
    l[i] <- .C("GetMPPModelName", as.integer(i-1),
               n=paste(rep(" ", .mpp.maxchar), collapse=""),
               PACKAGE="MarkedPointProcess")$n
  }
  return(l) 
}

splitmodel <- function(model) {
  if (missing(model) || (length(model)==0)) return(list(RF=list(),mpp=list()))
  op.list <- c("+","*")
  model.names <- get.mpp.names()
  all.names <- c(model.names, GetModelNames())

  trend <- model$trend
  mean <- model$mean
  model$trend <- model$mean <- NULL
  if (any(is.na(sapply(model, function(x) {
                        ifelse(is.null(x$m), pmatch(x, op.list),
                              pmatch(x$m, all.names)) }))))
    ## prevents that names are not given uniquely; e.g. model="n"
    ## could mean "nearest neighbour" or "nugget"
    stop(paste("operators not correct, or model names not unique within {",
               paste(all.names,collapse=", "),"}"))
  n <- sapply(model, function(x) ifelse(is.null(x$m), NA,
                                        pmatch(x$m, model.names)))
  index <- !is.na(n)
  index.op <- c(index[-1],FALSE) | c(FALSE,index[-length(index)])
  ## immer noch nicht 100% geprueft:
  RF <- model[!index.op]
  RF[index[!index.op]] <- "+"
  if ((length(RF)>0) && (RF[[1]]=="+"))  RF <- RF[-1]
  if ((length(RF)>0) && (RF[[length(RF)]]=="+")) RF <- RF[-length(RF)]
  if ((length(model)>1) &&
      any(unlist(model[index.op | c(FALSE,index[-length(index)])])!="+"))
    stop("marked point process models can be combined only additively") 
  mpp <- list()
  if (!(all(is.na(n)))) {
   n <- n[!is.na(n)] - 1
   np <-  integer(length(n))
    .C("GetNrMPPParameters", as.integer(n), as.integer(length(n)), np,
       PACKAGE="MarkedPointProcess", DUP=FALSE)
    model <- model[index]
    for (i in 1:length(model)) {
      if (length(model[[i]]$p)!=np[i]){
        print(c(length(model[[i]]$p),i,np[i]))
        stop("number of parameters incorrect")
      }
      mpp[[i]] <- list(model=model[[i]]$m, param=model[[i]]$p, mnr=n[i])
    }
  }
  if (!is.null(mean)) RF$mean <- mean
  if (!is.null(trend)) RF$trend <- trend  
  return(list(RF=RF,mpp=mpp))
}

simulateMPP <- function(coordmodel=c("given", "uniform", "Poisson"),
                         coord=NULL, npoints=NULL, lambda=NULL,
                         window=NULL, edgecorrection=0,
                         repetitions=1, coordrepet=1, model=NULL,
                         register=0, method=NULL)
{
  ## coord : created by coordmodel, if coordmodel>0 else coord must be given
  ##         n x 2 matrix
  ## npoints : how many points should be created at most, see coordmodel
  ## coordmodel: 0 : coord must be given
  ##             1 : npoints points, uniformly distributed in
  ##                 coordparameter[1,2] x coodparameter[3,4]
  ##             2 : Poisson distribution on
  ##                 coordparameter[1,2] x coodparameter[3,4] with
  ##                 lambda=coordparameter[6]; maximal npoints points
  ## window : 1,2 : xlim;  3,4:ylim
  ## repetitions    : number of realisations of marks with the same coord
  ## edgecorrection : interesting if nearest neighbour or additive boolean;
  ##                  Poisson distributed points are simulated outside area in
  ##                  frame of width coodparameter[5]
  ## coordrepet : number of realisations of the coord process;
  ##              is ignored if coordmodel=0
  ##
  ## model      : list of models a la InitSimulateRF, v2.0.0 onwards
  ##              new model names: nearest neighbour, random coin
  ##
  ## register,
  ## method     : see GaussRF()
  ##
  ## MPP-models
  ## ==========
  ## nearest neighbour: 1st. param. : factor with which the neirest neighbour
  ##                                  distance are multiplied
  ## random coin      : 1st param : 1 : disks,  2 : cone
  ##                    2nd param : radius of the disk or cone
  ##                    3rd param : height of the disk or cone

  PrintLevel <- RFparameters()$PrintLevel
  coordmodel <- pmatch(coordmodel, as.character(formals()$coordmodel)[-1])
  switch(coordmodel,
         { # given
           coordrepet <- 1
           if (is.null(npoints)) npoints <- nrow(coord)
           else stopifnot(npoints==nrow(coord))
           stopifnot(ncol(coord)==2)
           if (is.null(window)) window <- c(range(coord[,1]), range(coord[,2]))
           coord <- as.double(coord)
         },
         { # uniform
           stopifnot(!is.null(npoints), length(window)==4, is.null(coord))
           coord <- double(2 * npoints)
         },
         { # Poisson
           stopifnot(is.null(npoints), length(window)==4, is.null(coord))
           npoints <- 10 + qpois(0.9999999999, lambda = lambda *
                                 diff(window[c(1,2)]) * diff(window[c(3,4)]))
           ## good enough upper bound -- should never be reached
           coord <- double(2 * npoints)
         }
         ) 
  stopifnot(coordrepet>=1, repetitions>=1)
  if (edgecorrection>0 && is.null(lambda))
    lambda <- npoints / (diff(window[c(1,2)]) * diff(window[c(3,4)]))
  model <- splitmodel(model)
  RF <- model$RF
  mpp.models <- unlist(sapply(model$mpp, function(x) x$mnr))
  mpp.parameters <- unlist(sapply(model$mpp, function(x) x$p))
 
  error <- integer(1)

   data <- double(npoints * repetitions)
  result <- list();
  
  for (r in 1:coordrepet){
    result[[r]] <- list()
    act.npoints <- as.integer(npoints)
    ## must stay inside the loop
    ## since act.npoints is changed by GenerateMPPData
    
    .C("GenerateMPPData",
       coord,
       act.npoints,
       as.integer(coordmodel - 1),
       as.double(edgecorrection),
       as.double(window),
       as.double(lambda),
       data,
       as.integer(repetitions),
       as.integer(mpp.models),
       as.integer(length(mpp.models)),
       as.double(mpp.parameters),
       as.integer(length(mpp.parameters)),
       as.integer(PrintLevel),
       error, PACKAGE="MarkedPointProcess", DUP=FALSE)

    if (error!=0) stop(paste("Error in simulateMPP", error))
    ## IMPORTANT! if coordmodel is 2, then actual and physical act.npoints do
    ## not match; however matching is needed in the following steps

    result[[r]]$coord <- matrix(coord, ncol=2)[1:act.npoints,,drop=FALSE]

    result[[r]]$data <- matrix(data, ncol=repetitions)[1:act.npoints,,drop=FALSE]
    
    if (length(RF)>0) {
      rf <- GaussRF(x=result[[r]]$coord,
                    grid=FALSE,
                    model=RF,
                    n=repetitions,
                    method=method)
      if (is.null(rf)) stop("Error in random field simulation")
      result[[r]]$data <- result[[r]]$data + rf
    }
  }
 if (coordrepet==1) return(result[[1]]) else return(result)
}

         

rfm.test <- function(coord=NULL, data, normalize=TRUE,
                     MCrepetitions=99,  
                     MCmodel=list(model="exponential",
                       param=c(mean=0,variance=NA,nugget=0,scale=NA)),
                     method=NULL,
                     bin=c(-1,seq(0,1.2,l=15)),
                     MCregister=1,
                     n.hypo=1000,
                     pvalue=c(10, 5, 1),
                     tests="l1 & w3",
                     tests.lp=NULL, tests.weight=NULL,
                     Barnard=FALSE,
                     PrintLevel=RFparameters()$Print,
                     ...
                     )
{
  if (any(idx <- pvalue > 50)) {
    warning("old definition of the pvalue has beccome obsolete. 100 - pvalue is here used instead")
    pvalue <- 100 - pvalue
  }
  old.rf <- RFparameters()[c("Storing", "PrintLevel")]
  on.exit({if (!old.rf$Storing) DeleteRegister(MCregister);
           RFparameters(old.rf)})
  RFparameters(Storing=TRUE, PrintLevel=PrintLevel)
  
  stopifnot(n.hypo>1)
  if (PrepareModel(MCmodel)$anisotropy)
    stop("anisotropic models are not allowed (yet)")
  ## coord : coordinates
  ##         NULL  : full information is given by data==list(list(coord,data))
  ##         matrix: single data=marks set to analyse
  ##         list of matrices: several data sets, whose results are added 
  ## data  : vector or matrix or
  ##         list(list(coord,data)), but then length(data)==1, necessarily.
  ## normalize: should the data be normalizes first?
  ## MCrepet : number of simulations the MC-test is based on
  ## MCmodel : model is a list or double list, see CovarianceFct
  ##           [[ param : c(mean,variance,nugget,scale,further.parameters).
  ##              The ones that are NULL are estimated ]]
  ## sill  : if not NA, variance and nugget must be NA and the relation
  ##         variance+nugget=sill holds
  ## bin : bins for the "variograms"
  ##       excludes left margin and includes right margin of the bin
  ##
  ## Barnard : the permutation test by Barnard & Besag/Diggle

  ## if additive a table is returned with MCrepet + 1 columns
  ## else E/VAR/SD-test values normed to 100 (so as if MCrepet was 99) are
  ##      returned

  ## fitvario takes most of the time!

  distrNr <- as.integer(pmatch("Gauss", GetDistributionNames()) - 1) 
  ## maybe others will be allowed in future
  dummy <- .C("GetmppParameters", lnorms=integer(1), weights=integer(1),
                     tests=integer(1), mppmaxchar=integer(1), modelnr=integer(1),
                     PACKAGE="MarkedPointProcess", DUP=FALSE)
  .mpp.tests <- dummy$tests
  .mpp.weights <- dummy$weights
  .mpp.weightnames2 <- paste("w", 1:.mpp.weights, sep="")
  dummy <- as.matrix(expand.grid(.mpp.lpnames, .mpp.weightnames))
  .mpp.testnames <- c(paste(dummy[,1], " & ", dummy[,2], sep=""),
                      .mpp.extranames)
  dummy <- as.matrix(expand.grid(.mpp.lpnames, .mpp.weightnames2))
  .mpp.testnames2 <- c(paste(dummy[,1], " & ", dummy[,2], sep=""))
  
  if (tests=="all") users.tests <- rep(TRUE, .mpp.tests)
  else {
    users.tests <- rep(FALSE, .mpp.tests)
    if (!is.null(tests.lp) && !is.null(tests.weight)) {
      tw <- pmatch(tests.weight, .mpp.weightnames)
      tw[is.na(tw)] <- pmatch(tests.weight[is.na(tw)], .mpp.weightnames2)
      tests <- c(tests, pmatch(outer(pmatch(tests.lp, .mpp.lpnames),
                                     mpp.l.norms * (tw -1), "+")))
    }
    if (!is.null(tests)) {
      users.tests[pmatch(tests, .mpp.testnames)] <- TRUE
      users.tests[pmatch(tests, .mpp.testnames2)] <- TRUE
    }
    if (length(tests) > sum(users.tests)) {
      if (PrintLevel>0) {
        cat("Choose within:\n")
        cat(paste("'",.mpp.testnames,"'", sep=""), sep=", ")
        cat("\nor choose within:\n")
        cat(paste("'", c(.mpp.testnames, .mpp.extranames), "'",
                                 sep=""), sep=", ")
        cat("\n")
      }
      warning("some tests could not be matched or have been given twice")
    }
  }
  
  if (!is.null(coord)) {
    data <- list(list(coord=coord,data=data))
  } else {
    if (!is.null(data$data)) data <- list(data)
  }

  data <- lapply(data, function(x) list(coord=x$coord, data=as.matrix(x$data)))
  additive <- length(data)>1 || ncol(data[[1]]$data) > 1  
  lEbinM1 <-  as.integer(if (additive) MCrepetitions + 1 else 1)
  
  Etest <- integer(lEbinM1 * .mpp.tests)
  VARtest <- integer(lEbinM1 * .mpp.tests)
  SQtest <- integer(lEbinM1 * .mpp.tests) ###
  MAXtest <- integer(lEbinM1 * .mpp.nr.maxtests) ###
  error <- integer(1)

  norm <- function(x) {
    idx <- is.finite(x)
    if (any(idx)) x[idx] <- qnorm((rank(x[idx]) - 0.5) / sum(idx))
    return(x)
  }
  
  est <- list()
  for (r in 1:length(data)) {
    Data <- data[[r]]$data
    if (normalize) Data <- apply(Data, 2, norm)
    est <- list() ## only the last estimations are returned -- otherwise
    ##               it is too much
    
    for (i in 1:ncol(Data)) {
      idx <- !is.na(Data[, i])
      x <- data[[r]]$coord[idx, , drop=FALSE]
      d <- Data[idx, i]
      ld <- length(d)

      est[[i]] <- fitvario(x=x, y=NULL, z=NULL, T=NULL,
                           data=d, mle.methods="ml",
                           model=MCmodel, cross.methods=NULL, ...)$variogram$ml
      
      if (MCrepetitions>0) {     
        simu <- if (Barnard) {
          d[sample(ld, size=ld * MCrepetitions, replace=TRUE)]
        } else {
          GaussRF(x=x, grid=FALSE, model=est[[i]],
                  method=method, register=MCregister, n=MCrepetitions)
        }
        storage.mode(simu) <- "double"
       
        # if (PrintLevel>3) cat("Barnard-Besag-Diggle test on no correlation among marks (between marks and locations)!\n");

        
        .C("MCtest",  ## achtung! MCtest auch in tests/CHECK.R verwendet !
           as.integer(MCrepetitions),	
           as.double(x),
           as.double(d),
           as.integer(ld),
           as.integer(ncol(x)),
           simu,
           as.integer(PrintLevel),
           as.double(bin),     ## bins for the E & V functions
           as.integer(length(bin)-1),
#           as.double(Ebin),    ## bins for the ranks of the tests (in percent)
#           as.integer(lEbinM1),
           Etest,
           VARtest,
           SQtest,  ###
           .mpp.maxtests,
           .mpp.nr.maxtests,
           MAXtest,
           error,
           as.integer(additive), #as.integer(1),#additive=TRUE,
           ##                    28.8.04: rueckgeaendert; CHECK !!!!!
           as.integer(0), # do not copy simulation results & random choice
           ##               in mcf_internal
           PACKAGE="MarkedPointProcess", DUP=FALSE, NAOK=TRUE
           )
          if (error) stop(paste("error ",error))
      } # MCrep > 0
    } # ncol Data
  } ## 1..length(data)
                                     
  null.hypo <- null.sl <- reject.null <- NULL ## sl: significance level
  if (!additive) {
    null.sl <- list()
    reject.null <- list()
    if (!is.null(pvalue)) {
      stopifnot(length(pvalue) > 0)
      data.sl <- list(Etest, VARtest, SQtest)
     
      idx <- !is.na(data[[1]]$data[, 1])          
      x <- data[[r]]$coord[idx, , drop=FALSE]
      simu.data <- GaussRF(x=x, grid=FALSE, model=est[[1]],
                           method=method, n=n.hypo, register=MCregister)
                           
      null.hypo <-
        rfm.test(coord = x,
                 data=simu.data,
                 normalize=normalize,
                 MCrepetitions = MCrepetitions, MCmodel = MCmodel,
                 method=method,
                 bin = bin, 
                 MCregister = MCregister,
                 pvalue=NULL, ## otherwise risk of endless loop
                 tests="all",
                 PrintLevel=PrintLevel,...
                 )

      ## significance level of the null hypothesis for given pvalue     
      for (i in 1:length(data.sl)) { #E, Var, SQ: keep ordering in return below
        rate <- apply(null.hypo[[i]], 2, cumsum) #each column of null.hypo[[i]]
        ## contains the number of simulations where the null.sl of
        ## the simulated data set among the MC test simulations has been
        ## n, n=1,...,100 (if 100 MC test simulations have been performed)
        ## and the the matrix null.hypo[[i]] has 100 rows
        
        null.sl[[i]] <- matrix(nrow=length(pvalue), ncol=.mpp.tests)
        for (p in 1:length(pvalue)) {
          abspvalue <- (1 - pvalue[p] / 100) * n.hypo
          null.sl[[i]][p,] <- colSums(rate <= abspvalue) + 1 #
        }
   
        ## second: comparison with the estimated significance levels for the data
        reject.null[[i]] <- null.sl[[i]] <= rep(data.sl[[i]], eac=length(pvalue))
        
        dimnames(null.sl[[i]]) <- dimnames(reject.null[[i]]) <-
          list(pvalue, .mpp.testnames)
        reject.null[[i]] <- reject.null[[i]][, users.tests, drop=FALSE]
        null.sl[[i]] <- null.sl[[i]][, users.tests, drop=FALSE]

        if (any(null.sl[[i]] == 1))
          warning("strange result: null.sl equals 1 somewhere")
        if (any(null.sl[[i]] > n.hypo) && PrintLevel > 1)
          print("estimated position for p-value in null hypothesis exceeds position given by 'pvalue' by leading to imprecise results (",i,",",p,",",z,") -- increase MCrepetitions!\n")
      }
      names(reject.null) <- names(null.sl) <- c("E", "VAR", "SQ")
    }
  } else {
    if (!is.null(pvalue))
      warning("pvalue is not NULL, but unused (too many data sets) -- set pvalue=NULL?") 
  }
  Etest <- matrix(Etest, nrow=lEbinM1)
  VARtest <- matrix(VARtest, nrow=lEbinM1)
  SQtest <- matrix(SQtest, nrow=lEbinM1)
  dimnames(Etest) <- dimnames(VARtest) <- dimnames(SQtest) <-
    list(NULL, .mpp.testnames)
  
  return(list(E=Etest[, users.tests, drop=FALSE],     ## must be first 
              VAR=VARtest[, users.tests, drop=FALSE], ## must be second
              SQ=SQtest[, users.tests, drop=FALSE],   ## must be third
              ##  M = matrix(MAXtest,nrow=lEbinM1),
              ##  M not checked yet, so not returned yet
              reject.null=reject.null,
              est= est, ## only the results of the last element of the data list 
              ## MT= .mpp.maxtests, ## obsolete since this equals ncol(E)
              normalize=normalize, MCrepetitions=MCrepetitions,
              MCmodel=MCmodel,
              null.hypo=null.hypo,
              null.sl=null.sl,
              bin=bin))
}
  

######################################################################
######################################################################
## formerly: GetEfunction
mpp.characteristics <- function(...,
                                bin=NULL, rep=1, p=0.8, name="", normalize=TRUE,
                                show=FALSE, model=NULL, param=NULL,
                                summarize=TRUE,
                                PrintLevel=RFparameters()$Print,
                                dev=if (name=="") 2 else FALSE,
                                rdline=if (is.logical(dev)) NULL else readline,
                                staticchoice=FALSE){
  ## bin: c(-1,0,...,right margin of the last bin) for the E and VAR-functions
  ##      the programme needs "c(-1,0," as the first vector values!
  ##      (otherwise just nonsense is obtained as result for the tests)
  ## rep: data can be multivariate or repeated;
  ##      if they are both, the sequence is
  ##      (multivariate marks of repetition 1),  (multivariate marks of
  ##      repetition 2),...
  ##      data is a matrix(nrow=#individuals, ncol=(dim of multivariate marks)
  ##                                                                 * repetition
  ##      the program needs rep as an indicator, how ncol is split.
  ## p  : the outlier threshold for the robust estimations of the test statistics
  ## ...: a sequence of (matrix(ncol=2) of coordinates ,data), i.e.,
  ##       coord <- species1, marks <- species1, coord <- species2,
  ##                                                        marks <- species2,...
  ## name: if "" graphics are printed on screeen, otherwise stored in file named
  ##       name.*.*.*.ps
  ## normalize : transformation to marginal Gaussian data
  ## show: graphics are shown/printed or not
  ## model, param: a variogram model to compare with
  ## cex*: see ?par
  ##
  ## IMPORTANT NOTE: when calling this function,
  ##       either *all* the above parameters should be given explicitly
  ##       or the following data in "..." should be named, i.e,
  ##         GetEfunction(bin=seq(0,1,0.1), rep=1, p=0.8, name="",normalize=TRUE,
  ##                      show=FALSE, est=NULL, modelnr=0, coord1, data1)
  ##       or
  ##         GetEfunction(bin=seq(0,1,0.1),c1=coord1,d1=data1)
  ##       if a parameter is misspecified, the programme consequently interprets
  ##       as input data and strange errors will occur, that are sometimes
  ##       difficult to detect!
  ##
  ## Results:
  ## ========
  ## see MPPanalyse.cc for the specific ordering of the values of k_mm and the
  ##     mark variogram!

  dummy <- .C("GetmppParameters", lnorms=integer(1), weights=integer(1),
                     tests=integer(1), mppmaxchar=integer(1), modelnr=integer(1),
                     PACKAGE="MarkedPointProcess", DUP=FALSE)
  .mpp.tests <- dummy$tests

  args <- list(...)
  if (length(args)==1) {
    if (PrintLevel>1) print("length args=1 -- assuming a list!");
    args <- args[[1]]
  }
  if (PrintLevel>2)
    print(paste(c("Arguments:", names(args), "."), collapse=" "), quote=FALSE)
  
  n <- as.integer(length(args) / 2)  ## number of species
  nbins <- as.integer(length(bin)-1)
  dim <- as.integer(2);
  p <- as.double(p);
  bin <- as.double(bin)
  if ((col <- as.integer(ncol(as.matrix(args[[2]]))/rep))<1)
    stop("rep incorrect")

  rep <- as.integer(rep)
  col2 <- col *(col+1)/2
  colsum2 <- (n*col)*(n*col+1)/2


  for (i in ((1:n)*2)) {
    args[[i]] <- as.matrix(args[[i]])
    if ( !(all(apply(args[[i]], 2, var, na.rm=TRUE) != 0)) ) {
      stop("data set seems to be a trivial one")
    }
  }

  if (normalize) {
    norm <- function(x) {
      idx <- is.finite(x)
      if (any(idx)) x[idx] <- qnorm((rank(x[idx]) - 0.5) / sum(idx))
      return(x)
    }
    for (i in ((1:n)*2)) {
      ## originally not only for each species and each variable, but also
      ## also for each repetition separately:
      ##
      args[[i]]  <- apply(args[[i]], 2, norm)
      ##
      ## now repetetions are put together
      ## note: apply(...,1,norm) transforms matrix at the same time!
    }
  }
    
  if (is.null(bin)) {
    ## choose maximal distance as max bin end,
    ## number of bins = min <- species(#indiviuals) / 5
    stop("null bin not programmed yet")
  }
  ## investigate if all the data have the same number of columns; if not, stop
  ## and if divisible by rep
  E <- double(rep * n^2 * col * nbins)
  ETest <- double(rep * n^2 * col * .mpp.tests)
  Ebin <- integer(n^2 * col * nbins)
  VAR <- double(rep * n^2 * col2 * nbins)
  VARTest <- double(rep * n^2 * col2 * .mpp.tests)
  ##  MAXTest <- double(rep * n^2 * col2 * .mpp.nr.maxtests)
  SQ  <- double(rep * n^2 * col2 * nbins)
  SQTest <- double(rep * n^2 * col2 * .mpp.tests)
  VARbin <- integer(n^2 * col2 * nbins)
  KMM <- double( rep * colsum2 *  nbins)
  KMMbin <- integer(colsum2 * nbins)
  GAM <- double(rep * colsum2 *  nbins)
  GAMbin <- integer(colsum2 * nbins)
  error <- integer(1);
  PrintLevel <- as.integer(PrintLevel)

  
  block <- ""
  for (i in (1:n)) {
    block <- paste(block,"as.double(args[[",2*i-1,
          "]]),as.double(args[[",2*i,
          "]]),as.integer(nrow(args[[",2*i-1,"]])),",sep="");
  }
  
  text<-paste('.C("mcf",',
              "E,ETest,Ebin,VAR,VARTest,SQ,SQTest,VARbin,KMM,KMMbin,GAM,GAMbin,",
              "error,PrintLevel,p,bin,nbins,dim,n,col,as.integer(staticchoice),",
              "rep,", block,"PACKAGE='MarkedPointProcess',DUP=FALSE,NAOK=TRUE)",
              sep="")
  # print(text)
  eval(parse(text=text))

  if (summarize && (rep>1)) {
    E <- rowMeans(matrix(E, ncol=rep))
    ETest <- rowMeans(matrix(ETest, ncol=rep))
    VAR <- rowMeans(matrix(VAR, ncol=rep))
    VARTest <- rowMeans(matrix(VARTest, ncol=rep))
    SQ <- rowMeans(matrix(SQ, ncol=rep))
    SQTest <- rowMeans(matrix(SQTest, ncol=rep))
    KMM <- rowMeans(matrix(KMM, ncol=rep))
    GAM <- rowMeans(matrix(GAM, ncol=rep))
    ## MAXTest <- rowMeans(matrix(,ncol=rep))
    rep <- 1
  }

  E <- matrix(E, nrow=nbins)
  ETest <- matrix(ETest, nrow=.mpp.tests)
  VAR <- matrix(VAR, nrow=nbins)
  VARTest <- matrix(VARTest, nrow=.mpp.tests)
  SQ <- matrix(SQ, nrow=nbins)
  SQTest <- matrix(SQTest, nrow=.mpp.tests)
  KMM <- matrix(KMM, nrow=nbins)
  GAM <- matrix(GAM, nrow=nbins)
  ## MAXTest <-
  midbin <- c(0,(bin[c(-1,-length(bin))]+bin[c(-1,-2)])/2)
  if (show) {
    k <- 1;
    l <- 1
    rangebin <- c(0,bin[length(bin)])
    ps <- paste(name, "Ebin.ps", sep=".")
    Dev(TRUE, dev, ps=ps)
    plot(midbin, Ebin, xlab="distance r", ylab="n")
    Dev(FALSE, dev)
    if (!is.null(rdline)) rdline(ps)
    ps <- paste(name, "VARbin.ps", sep=".")
    Dev(TRUE, dev, ps=ps)
    plot(midbin, VARbin, xlab="distance r", ylab="n")
    Dev(FALSE, dev)    
    for (i in 1:n) { #species i
      segmentE <- (i - 1) * (col * n + 1) 
      segmentV <- (i - 1) * (col2 * n + 1)
      for (j in 1:n) { # species j
        for (cc in 1:col) {
          ## ** E of mark cc of species i given species j **
          ylab <- "E";
          if (n!=1) ylab <- paste(ylab," <- {",i,"|",j,"}",sep="")
          ylab <- paste(ylab,"(",sep="")
          if (col!=1) ylab <- paste(ylab,"m <- ",cc,";",sep="")
          ylab <- paste(ylab,"r)",sep="")
          for (r in ((0:(rep-1)) * col * n * n)){
            if (!is.null(rdline)) rdline(ps)
            ps <- paste(name, "E", i, j, cc, "ps", sep=".")
            Dev(TRUE, dev, ps=ps)
            plot(midbin, E[, k + r], xlab="distance r", ylab=ylab)
            lines(rangebin, c(E[1, r + segmentE + cc],
                              E[1, r + segmentE + cc]), lty=3)
            Dev(FALSE, dev)
          }
          k <- k+1;          
        }
        ll <- l
        for (cc in 1:col) {
          for (cc2 in cc:col) {
            ## ** Cov(cc, cc2) of species j given species i **
            ## i varies slowest, then j, then mark cc, then mark cc2 (restricted)
            ylab <- ifelse(cc==cc2, "V", "V");
            if (n!=1) ylab <- paste(ylab," <- {",i,"|",j,"}",sep="")
            ylab <- paste(ylab,"(",sep="")
            if (col!=1) ylab <- paste(ylab,"m_",cc,"m_",cc2,";",sep="")
            ylab <- paste(ylab,"r)",sep="")
            
            for (r in ((0:(rep-1)) * col2 * n * n)){
             if (!is.null(rdline)) rdline(ps)
             ps <- paste(name, "V", i, j, cc, "ps", sep=".")
             Dev(TRUE, dev, ps=ps)
              plot(midbin,VAR[, l+r],xlab="distance r",ylab=ylab)
              lines(rangebin,c(VAR[1, r+segmentV + l - ll + 1],
                               VAR[1, r+segmentV + l - ll + 1]), lty=3)
              Dev(FALSE, dev)
            }
            l <- l + 1;
          }
        }
        l <- ll
        for (cc in 1:col) {
          for (cc2 in cc:col) {
            ylab <- ifelse(cc==cc2, "SD", "Cor");
            if (n!=1) ylab <- paste(ylab, " <- {", i, "|", j, "}", sep="")
            ylab <- paste(ylab,"(",sep="")
            if (col!=1) ylab <- paste(ylab, "m_", cc, "m_", cc2, ";", sep="")
            ylab <- paste(ylab,"r)", sep="")
            
            for (r in ((0:(rep-1)) * col2 * n * n)){
              if (!is.null(rdline)) rdline(ps)
              ps <- paste(name, "SQ", i, j, cc, "ps", sep=".")
              Dev(TRUE, dev, ps=ps)
              plot(midbin, SQ[, l + r], xlab="distance r", ylab=ylab)
              lines(rangebin,c(SQ[1, r + segmentV + l - ll + 1],
                               SQ[1, r + segmentV + l - ll + 1]),
                    lty=3)
              Dev(FALSE, dev)
            }
            l <- l + 1;
          }
        }
      } # j
    } # i
    
    if (!is.null(rdline)) rdline(ps)
    ps <- paste(name, "KMMbin.ps", sep=".")
    Dev(TRUE, dev, ps=ps)
    plot(midbin, KMMbin, xlab="distance r", ylab="n")
    Dev(FALSE, dev)    
    k <- 1;
    for (i in 0:(n*col-1)) {
      ## (cross) k_mm and (cross) mark variogram
      ## see MPPanalyse.cc for the specific ordering of these values!
      spec1 <- as.integer(i/col) +1;
      col1  <- i - (spec1-1)*col;
      for (j in i:(n*col-1)) {
        spec2 <- as.integer(j/col) +1;
        col2 <- j - (spec2-1)*col;

        ylab <- ifelse(n==1, "", paste("_{",spec1,",",spec2,"}",sep=""))
        ylab <- paste(ylab,"(",sep="")
        if (col!=1) ylab <- paste(ylab,"m_",col1,"m_",col2,";",sep="")
        ylab <- paste(ylab,"r)",sep="")
        
        for (r in ((0:(rep-1))*colsum2)){
          if (!is.null(rdline)) rdline(ps)
          ps <- paste(name,"K",spec1,spec2,col1,col2,"ps",sep=".")
          Dev(TRUE, dev, ps=ps)
          plot(midbin,KMM[,k+r],xlab="distance r",ylab=paste("KMM",ylab,sep=""))
          Dev(FALSE, dev)
        }        
        for (r in ((0:(rep-1))*colsum2)){
          if (!is.null(rdline)) rdline(ps)
          ps <- paste(name,"G",spec1,spec2,col1,col2,"ps",sep=".")
          Dev(TRUE, dev, ps=ps)
          plot(midbin, GAM[,k+r], xlab="distance r",
               ylab=paste("gamma", ylab, sep=""))
          if (!is.null(model)) {
            xx <- seq(max(midbin)/100000,max(midbin),l=length(midbin)*4);
            lines(xx, Variogram(x=xx, model=model, param=param))
           }
        Dev(FALSE, dev)
        }        
        for (r in ((0:(rep-1))*colsum2)){
          if (!is.null(rdline)) rdline(ps)
          ps <- paste(name, "GvsK", spec1, spec2, col1, col2,"ps",sep=".")
          Dev(TRUE, dev, ps=ps)
          plot(KMM[,k+r],GAM[,k+r],xlab=paste("KMM",ylab,sep=""),
              ylab=paste("gamma",ylab,sep=""))
          Dev(FALSE, dev)
         }        
        k <- k+1
      }
    }      
  }
  ret <- list(E=E,ETest=ETest,VAR=VAR,VARTest=VARTest,
              SQ=SQ, SQTest=SQTest,
              KMM=KMM,GAM=GAM,
              Ebin=Ebin,VARbin=VARbin,KMMbin=KMMbin,GAMbin=GAMbin,
              midbin=midbin,
              call=match.call())
  if (show) invisible(ret) else return(ret)
}

