
##########################################################################
## pictures.R
##########################################################################

srd.jrssb <- function(input=NULL, repet=500, dev=2, PrintLevel=2, readlines=TRUE,
                      ps.path="ps/", data.save.path="data/", simu.path="simu/",
                      tex.path="tex/", biondi.etal=NULL,
                      final=TRUE
                      ) {
  on.exit(traceback())
  debug <- FALSE
  texfilename <- "srd.jrssb.04"
  systempath <- system.file(package='MarkedPointProcess')
  ## systempath <- "/home/schlather/R/MPP/MarkedPointProcess/inst"
  if (final) {
    optim.control <- NULL
    data.n.hypo <- 1000 ## number of simulations under the null hypothesis to 
    ##                get the significance level under the null hypothesis
    ##                for the given pvalue
    simu.gauss.variance <- c(0.000001, (1:10)/10)
    simu.R.vario <- c(0.0001, 0.01, 0.05, 0.1, 0.2)
    simu.radii <- c(0.05, 0.1, 0.15, 0.2, 0.25) 
    simu.R.nn.vario <- c(0.0001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.8)
    simu.R.nn.point.setting <-
      list(c(i=12.5, x=1, y=1),
           c(i=25, x=sqrt(2), y=sqrt(2)),
           c(i=25, x=1, y=1),
           c(i=50, x=1, y=1),
           c(i=50, x=sqrt(2), y=sqrt(2)),
           c(i=50, x=2, y=2),
           c(i=100, x=sqrt(2), y=sqrt(2)),
           c(i=100, x=2, y=2),
           c(i=100, x=2*sqrt(2), y=2*sqrt(2)),
           c(i=200, x=2, y=2),
           c(i=200, x=2*sqrt(2), y=2*sqrt(2)),
           c(i=200, x=2*2, y=2*2))
    simu.individ.list <- c(50, 100, 200)
  } else {
    #if (repet>3) {
    #  warning("repet>3 replaced by repet=3")
    #  repet <- 3
    #}
    optim.control <- list(factr=1e14)
    data.n.hypo <- 2 # 5
    simu.gauss.variance <- c(0.000001, (1:2)/2)
    simu.R.vario <- c(0.0001, 0.05) # do not change
    simu.radii <- c(0.05, 0.1)      # do not change
    simu.R.nn.vario <- c(0.0001, 0.05)
    simu.R.nn.point.setting <-
      list(c(i=12.5, x=1, y=1),
           c(i=25, x=sqrt(2), y=sqrt(2)),
           c(i=25, x=1, y=1),
           c(i=50, x=1, y=1),
           c(i=50, x=sqrt(2), y=sqrt(2)),
           c(i=50, x=2, y=2),
           c(i=100, x=sqrt(2), y=sqrt(2)),
           c(i=100, x=2, y=2),
           c(i=100, x=2*sqrt(2), y=2*sqrt(2)))
    simu.individ.list <- c(100)
  }
  
  data(BITOEK)
  RFparameters(Print=PrintLevel, pch=if (PrintLevel<3) "" else "*",
               TBMCE.force=TRUE)

  steigerwald.showEfct <- TRUE
  steigerwald.normalize <- TRUE
  coulissenhieb.showEfct <- TRUE
  coulissenhieb.normalize <- TRUE
  biondi.etal.showEfct <- TRUE
  biondi.etal.normalize <- TRUE
  data <-
    list(list(set="steigerwald", name="Steigerwald", pfactor=5),
         list(set="coulissenhieb", name="Coulissenhieb", pfactor=5),
         list(set="biondi.etal", name="Gus Pearson", pfactor=0.01)
         )
  data.fcts <- c("E", "VAR", "GAM", "Ebin")
  data.fctnames <-  c("E", "V", "G", "n")
  data.xlim <- c(0, 40) #  c(0, 50)
  data.ylims <- cbind(c(-0.4, 0.1), 
                 c(0.6, 1.1), 
                 c(0, 1.1), 
                 c(0, 2500))
  data.xlab <- "distance"
  data.ylabs <- c("E", "V", expression(gamma), "n")
  data.MCmodel <- list(model="whittle", param=c(0, NA, NA, NA, NA))
  data.bin <- c(-1, seq(0, 50, 2))
  data.conf.repet <- 100
  data.confid.p <- 0.05
  
  theo.pch <- c(0, 4, 1)
  theo.lwd <- 3
  theo.len <- 100
  theo.ltyE <- 2
  theo.ltyV <- 1
  theo.ltyG <- 3
  theo.pchE <- 0
  theo.pchV <- 1
  theo.pchG <- 2

  simu.fig1.norm.select <- c("max", "l1")
  simu.fig1.norm.select <- c("max", "l1")
  simu.fig1.weight.select <- c("w1", "w2", "w3", "w4", "w5", "w6", "w7")
  simu.fig1.weight.select <- c("w3", "w5", "w7")
  simu.fig1.individ <- 100
  simu.rlist <- list(list(radius=0.05, scale=0.05, gauss.variance=0, simu=FALSE),
                     list(radius=0.1, scale=0.05, gauss.variance=0.5, simu=FALSE)
                     ) ## the values used in Fig 1, except for gauss.variance
  ##                 where standard is used, see below
  simu.standard.MC <- simu.standard.cov <- simu.standard.estim <- "exponential"
  simu.models <- list(add=list(file=paste("additive", repet, sep="-"),
                        title="Random coin model", model="random coin",
                        intrinsic.var=function(lbRd) lbRd),
                      nn=list(file=paste("nn", repet, sep="-"),
                        title="Nearest neighbour", model="nearest neighbour"),
                      var=list(file=paste("variance", repet, sep="-"),
                        title="Variance by random coins", model="variance",
                        intrinsic.var=function(lbRd) 1 + lbRd * (3 + lbRd)),
                      )
  assign(env=NULL, "simu.completename",
    function(basicname, variance) paste(basicname, variance * 100, sep="."))
  assign(env=NULL, "simu.addname",
    function(basic.name, estimate, simulate, MCmodel, gauss.model,
             radius, scale, individ)
      paste(basic.name, estimate, simulate, MCmodel, gauss.model, 
            radius*100, scale *100, individ, sep="."))
  assign(env=NULL, "simu.nnname",
         function(scale, x, y, individ)
         paste(simu.models$nn$file, scale * 100, round(x*y), individ, sep="."))
  simu.name.EVS <- c("E", "V", "S")
  simu.MCrep <- 99  ## function not tested for simu.MCrep != 99 !!

  if (is.numeric(dev)) {
    fig1.cex <- 1
    parlist <- list(cex=1, cex.axis=1, cex.lab=1, cex.main=1, lwd=3, 
                    cex.sub=1, mar=c(2.2, 2.6, 0.4, 0.4), pch=16, las=0)
  } else {
    fig1.cex <- 0.9
    parlist <- list(cex=1, cex.axis=1.5, cex.lab=1.5, cex.main=1.5, lwd=2, 
                    cex.sub=1.5, mar=c(2.2, 2.6, 0.4, 0.4), pch=16, las=0)
  }

  for (n in c("ps.path", "tex.path", "data.save.path", "simu.path")) {
    name <- get(n)
    if (PrintLevel>5) cat("path", name, "\n")
    if (file.exists(name)) stopifnot(file.info(name)$isdir)
    else stopifnot(dir.create(name))
    if (substr(name, nchar(name), nchar(name)) != "/")
      assign(n, paste(name,  "/", sep=""))
  }
  
  assign(env=NULL, "norm", function(x) {qnorm((rank(x)-0.5)/length(x))})

  assign(env=NULL, "rl", if (readlines)
         function(x) if (x=="") readline("press return")
         else readline(paste(x, ": press return"))
  else function(x) { cat(x); sleep.milli(1000); cat("\n")})
        
  .dummy <- .C("GetmppParameters", lnorms=integer(1), weights=integer(1),
               tests=integer(1), mppmaxchar=integer(1), modelnr=integer(1),
               debug=integer(1),
               PACKAGE="MarkedPointProcess", DUP=FALSE)
  .mpp.tests <- .dummy$tests
  .mpp.weights <- .dummy$weights
  .mpp.weightnames2 <- paste("w", 1:.mpp.weights, sep="")

  .dummy <- as.matrix(expand.grid(.mpp.lpnames, .mpp.weightnames))
  .mpp.testnames <- c(paste(.dummy[,1], " & ", .dummy[,2], sep=""),
                      .mpp.extranames)
  .dummy <- as.matrix(expand.grid(.mpp.lpnames, .mpp.weightnames2))
  .mpp.testnames2 <- c(paste(.dummy[,1], " & ", .dummy[,2], sep=""))
  .mpp.l.norms <- length(.mpp.lpnames)
  stopifnot(100 %in% simu.individ.list)

######################################################################
######################################################################

  
  simulate.rfmtest <- function(coord=NULL, 
                               npoints=NULL, lambda=NULL, 
                               coordmodel="uniform", window, 
                               repetitions=1, edgecorrection=0.0, coordrepet=1, 
                               model, register=0,
                               method=(if (is.null(npoints) ||  npoints>500)
                                       NULL else "direct"), 
                               gauss.variance=simu.gauss.variance, 
                               
                               normalize=TRUE, 
                               MCrepetitions=simu.MCrep, 
                               MCmodel= list(model=simu.standard.MC, 
                                 param=c(mean=0, variance=NA, nugget=0,
                                   scale=NA)), 
                               sill=NA, 
                               bin=c(-1, seq(0, 0.7, l=15)), 
                               Ebin=seq(0, 1, 0.01),           
                               use.naturalscaling=TRUE, 
                               Barnard=FALSE, 
                               MCregister=1, 

                               saving=NULL, path="./") {
    ## coord, ..., method, see simlate.mpp
    ## gauss.variance : part of the variance belonging to the Gaussian random 
    ##                  field (alpha in the section 4.1 of SRD (2004))
    ##
    ## normalize, ..., MCregister, see rfm.test
    ##
    ## saving : null : result is returned (not used in my studies)
    ##          true : result is saved in Result.*.*.* (not used in my studies)
    ##          string : saved under given name
    ##                   filenames should be chosen such that important
    ##                   running variables are reflected and the
    ##                   numerical values have been multiplied with 100
    ##                   see print.mpp.simu below and the worksheet files in
    ##                   this folder
    ##       if true or string the last part of the name is the respective value
    ##       of alpha (see gauss.variance above)
    ## path : for saving results, used in combination with saving
    ## output : main effect is the side effect of creating a file with the
    ##          simulated data

    if (PrintLevel>6)
      cat(npoints, lambda, coordmodel, window, coordrepet, MCmodel$model, "\n")
    stopifnot(is.function(model))
    
    if (is.character(saving)) {
      basicfn <- saving
      saving <- TRUE
    } else if (is.logical(saving))
      basicfn <- paste("mpp.result", npoints, coordmodel, 
                       coordparameter[length(coordparameter)], 
                       paste(model, collapse="_"), MCrepetitions,
                       if (is.function(MCmodel)) MCmodel(0)$model else MCmodel, 
                       variance*100, sep=".")
    else if (is.null(saving)) {
      basicfn <- ""
      saving <- FALSE
    } else stop("value of saving not allowed")
    basicfn <- paste(path, basicfn, sep="")
  
    for (variance in gauss.variance) {
      filename <- simu.completename(basicfn, variance)
      if (PrintLevel>7) cat(variance, model(variance)[[1]]$model, "\n")
      ## same simulation again?
      seed.file <- paste(filename, "seed", sep=".")
      
      if (FileExists(seed.file)) {
        if (PrintLevel>4) cat("loading", seed.file, "\n")
        load(seed.file)
        assign(".Random.seed", seed.save, envir=.GlobalEnv)
      } else {
        runif(1)
        seed.save <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
        if (saving) save(seed.save, file=seed.file)
        LockRemove(seed.file);
      }
      if (saving && FileExists(filename)) next
      
      data <- simulate.mpp(coord=coord, npoints=npoints, lambda=lambda,
                           coordmodel=coordmodel, window=window,
                           repetitions=repetitions,edgecorrection=edgecorrection,
                           coordrepet=coordrepet, model=model(variance),
                           register=register, method=method)
      result <- rfm.test(data=data, normalize=normalize, 
                         MCrepetitions = MCrepetitions,
                         MCmodel = (if (is.function(MCmodel)) MCmodel(variance)
                                    else MCmodel),
                         bin = bin, Ebin = Ebin, 
                         MCregister = MCregister, method=method, 
                         pvalue=NULL, tests="all", Barnard=Barnard,
                         sill = sill, use.naturalscaling = use.naturalscaling,
                         optim.control=optim.control)
      if (!debug) result$null.hypo <- NULL
      
      if (saving) {
        if (PrintLevel>3) cat("saving", filename, "\n")
        save(result, file=filename)
        LockRemove(filename);
      } else {
        if (length(gauss.variance)>1)
        warning("result only for the first entry of gauss.variance")
        return(result)
      }
      rm("data", "result")
    }
    return(NULL)
  } ## simulate.rfmtest
  
 
  plot.fig1 <- function(name, npoints,
                        gauss.variance=simu.gauss.variance, 
                        Title="", 
                        pvalue=c(90, 95), ## werden hier neu berechnet, und
                        ## nicht von rfm.test uebernommen
                        ps=NULL, 
                        absolute=TRUE, 
                        # PrintLevel=RFparameters()$Print, 
                        select.lp = "all", 
                        select.weight = "all", 
                        select.extra = FALSE, 
                        return.rate = FALSE, 
                        dev = if (is.null(ps)) 2 else FALSE, 
                        print.title = !is.logical(dev), 
                        height=5, width=5, 
                        coloured = TRUE, 
                        print.legend = TRUE, 
                        type="l", cex=1, 
                        y.pos.leg = 0.65,
                        pch.basic=1:.mpp.l.norms,
                        lty.basic=1:.mpp.weights,
                        lty.abline=4, col.abline=1, 
                        ...) {

    ## PURPOSE : printing of the summary results of a large simulation study
    ##           the variables given in ... are the variables running in
    ##           the file name
    ##           The algorithm is based on
    ##             * the assumption that the teststatistic is based on 99 MC
    ##               simulation
    ##             * Ebin=1:100
    ##
    ## plots parts of Figure 1 in SRD 04
    ##
    ## this function is still in a raw version!
    ## 
    ## name  : name of file
    ## npoints    : numbers of spatial points in each simulation
    ##                  typically 50, 100, 200 [vector]
    ## gauss.variance : results based on simulations where the parts [in %] of 
    ##                  the total variance caused by the random fields equals
    ##                  gauss.variance [vector]
    ## Title : string for the title of the plots
    ## pvalue: nominal p-values, i.e. the minimal position number of the
    ##         test statistic (ETest, VARTest, ...) to be counted as success 
    ## ps    : base name of the postscript files; if is.null then printed on the
    ##         screen
    ## absolute : if TRUE the raw test statistics are returned
    ##            if FALSE the test statistics are corrected taking the
    ##            estimations for gauss.variance=1 into account
    ##            (the usual correction introduced since parameters are
    ##            estimated for the null hypothesis -- here in a simplified
    ##            form (assuming that the errors introduced for different
    ##            values of alpha (gauss.variance) are about the same);
    ##            otherwise it would take too much time)
    ## PrintLevel : see RFparameters in RandomFields
    ## select.lp, select.weight : see tests.lp and tests.weight in rfm.test
    ##                            (except that the three additional test
    ##                             are not considered)
    ## select.extra : if TRUE than the three additional tests are included
    ## return.rate : if TRUE only the rates in the computation of the
    ##               correction for !absolut are returned (only for
    ##               debugging/developping)
    ## dev   : see Dev in RandomFields
    ## print.title, height, width, coloured, print.legend, cex, 
    ##       type : graphical parameters for the plots
    ## y.pos.leg : y position of the legend
    ## lty.basic : line type for the graphs 
    ## lty.abline, col.abline : line type and colour for the 5% level line
    simu.fig1.ps.special <- c("E", "V", "S", "M")
    users.tests <- rep(FALSE, .mpp.tests)
    if (length(select.lp)==1 && select.lp=="all") {
      sl <- 1:.mpp.l.norms
      select.lp <-rep(TRUE, .mpp.l.norms)
    } else {
      sl <- pmatch(select.lp, .mpp.lpnames)
      stopifnot(all(!is.na(select.lp)))
      select.lp <- rep(FALSE, .mpp.l.norms)
      select.lp[sl] <- TRUE
    }
    if (length(select.weight)==1 && select.weight=="all") {
      sw <- 1:.mpp.weights
      select.weight <- rep(TRUE, length(.mpp.weightnames))
    } else {
      sw <- pmatch(select.weight, .mpp.weightnames)
      sw[is.na(sw)] <- pmatch(select.weight[is.na(sw)], .mpp.weightnames2)
      stopifnot(all(!is.na(select.weight)))
      select.weight <- rep(FALSE, length(.mpp.weightnames))
      select.weight[sw] <- TRUE
    }
    
    select <- rep(FALSE, .mpp.tests)
    select[outer(sl, (sw - 1) * .mpp.l.norms, "+")] <- TRUE
    if (select.extra)
      select[.mpp.tests :(.mpp.tests-length(.mpp.extranames)+1)] <- TRUE
    select.extra <- rep(select.extra, length(.mpp.extranames))
    select.weight <- c(select.weight, select.extra)
    select.lp <- c(select.lp, select.extra)
    
    if (coloured)
      col.leg <- c(rainbow(.mpp.l.norms), rep("black", length(.mpp.extranames)))
    else {
      col.leg <- rep(NA, .mpp.l.norms + length(.mpp.extranames))
      col.leg[select.lp] <- grey(seq(0, 0.8, length=sum(select.lp)))
    }
    col <- c(rep(col.leg[1:.mpp.l.norms], .mpp.weights), 
             col.leg[-1:-.mpp.l.norms])[select]
    col.leg <- col.leg[select.lp]

    pch.leg <- c(pch.basic, .mpp.l.norms + 1:length(.mpp.extranames))
    pch <- c(rep(pch.leg[1:.mpp.l.norms], times=.mpp.weights), 
             pch.leg[-1:-.mpp.l.norms])[select]
   
    lty.leg <- c(lty.basic, rep(9, length(.mpp.extranames)))
    lty <- c(rep(lty.leg[1:.mpp.weights], each=.mpp.l.norms), 
             lty.leg[-1:-.mpp.weights])[select]
    
    Numbers <- c(rep(.mpp.tests, 3), .mpp.nr.maxtests) ## E1, V2, S3, MAX4
    Numbers <- rep(.mpp.tests, 3) ## E1, V2, S3 -- ONLY!
    
    args <- list(...);
    args <- args[sapply(args, function(y) !is.null(y))] ## remove NULLs
    bin <- c(-1, seq(0, 0.7, l=15))
    xlim <- c(0, 1)
    ylim <- c(0, 1)
    position <- list();
    factor <- list();
    for (i in 1:length(Numbers)) {
      position[[i]] <- matrix(nrow=length(pvalue), ncol=Numbers[i])
      factor[[i]] <- matrix(nrow=length(pvalue), ncol=Numbers[i])
    }
    index <- (1:length(gauss.variance))[gauss.variance == 1.0]
    if (absolute) {
      for (i in 1:length(Numbers)) {
        for (p in 1:length(pvalue)) {
          position[[i]][p, ] <- pvalue[p]
          factor[[i]][p, ] <- 1.0
        }
      }
    } else {
      if (length(index)!=1) stop("gauss.variance[]=1.0 not given or given twice")
      gauss.variance <- gauss.variance[-index]
    }
    
    without1 <- 1:length(gauss.variance)
    if (absolute && (length(index)>=1) && (length(gauss.variance)>1))
      {without1 <- without1[-index] }
    
    Xtotal <- list();
    for (i in 1:length(Numbers)) {
      Xtotal[[i]] <- list(length(pvalue))
      for (p in 1:length(pvalue)) {
        Xtotal[[i]][[p]] <-
          matrix(0, nrow=length(gauss.variance), ncol=Numbers[i])
      }
    }
    ENVIR <- environment()
    
    print.data <- function() {
      X <- list(length(Numbers)) ## E1, V2, S3, 
      for (i in 1:length(Numbers)) {
        X[[i]] <- list(length(pvalue))
        for (p in 1:length(pvalue)) {
          X[[i]][[p]] <- matrix(0, nrow=length(gauss.variance), ncol=Numbers[i])
        }
      }      
      
      if (!absolute) {
        name1 <- paste(basicname, 100, sep="")
        load(name1)
        for (i in 1:length(Numbers)) {
          rate <- apply(result[[i]], 2, cumsum)
          rate <- 1-rate/rate[nrow(rate), 1]           
          for (p in 1:length(pvalue)) {
            position[[i]][p, ] <- apply(rate>1-pvalue[p]/100, 2, sum)
            for (z in 1:Numbers[i]) {         
              if (position[[i]][p, z]>0)
                factor[[i]][p, z] <-
                  (rate[1+position[[i]][p, z], z]-(1-pvalue[p]/100))/
                    (rate[1+position[[i]][p, z], z]-rate[position[[i]][p, z], z])
              else {
                if (PrintLevel>2)
                  cat("position is zero at (", i, ", ", p, ", ", z, ")\n")
                stopifnot(is.na(factor[[i]][p, z]))
              }
            }
          }           
        }
      }
      number.files <- 0;
      for (v in 1:length(gauss.variance)) {
        name1 <- paste(basicname, gauss.variance[v]*100, sep="")
                                        #cat(v, ": ", name1, "\n")
        if (file.exists(name1)) {
          load(name1)
          if (gauss.variance[v]<1) {number.files <- number.files + 1;}
          for (i in 1:length(Numbers)) {
            rate <- apply(result[[i]], 2, cumsum)
            rate <- 1-rate/rate[nrow(rate), 1]
            if (absolute) {
              for (p in 1:length(pvalue)) {
                X[[i]][[p]][v, ] <- rate[pvalue[p], ] }
            } else {
              for (p in 1:length(pvalue)) {
                for (z in 1:Numbers[i]) {
                  if (position[[i]][p, z]>0)
                    X[[i]][[p]][v, z] <-
                      factor[[i]][p, z] * rate[position[[i]][p, z], z] +
                        (1-factor[[i]][p, z]) * rate[1+position[[i]][p, z], z]
                  else  {
                    if (PrintLevel>2)
                      cat("position is zero at (", i, ", ", p, ", ", z, ")\n")
                    stopifnot(is.na(factor[[i]][p, z]))
                  }
                }
              }
            }
          }
        } else {
          warning(paste("file", name1, "does not exist!"));
          for (i in 1:length(Numbers)) {
            for (p in 1:length(pvalue)) X[[i]][[p]][v, ] <- 0 
          }
        }
      } ## gausvariance
      if (return.rate) return(X)
      
      for (i in 1:length(Numbers)) {
        for (p in 1:length(pvalue)) {
          Xtotal[[i]][[p]] <-  Xtotal[[i]][[p]] + X[[i]][[p]]
          assign("Xtotal", Xtotal, envir=ENVIR)
          
          if (PrintLevel>5) print(c(i, p))
          filename <- paste(psbasicname, simu.fig1.ps.special[i],
                            pvalue[p], sep=".")
          Dev(TRUE, dev, height=height, width=width, ps=filename)
          percent <- 1
          percent <- 100
          
          matplot(gauss.variance, X[[i]][[p]][, select, drop=FALSE] * percent, 
                  ylim=c(0, percent), col=col, pch=pch, lty=lty, lwd=3, 
                  ylab=paste("rejection rate [%]"), #simu.fig1.ps.special[i]
                  xlab=expression(alpha), type=type, cex=cex)
          
          if (print.title) title(paste(Titel, "  p-value=", pvalue[p], sep=""))
          abline((1 - pvalue[p]/100) * percent, 0, col=col.abline,
                 lty=lty.abline, lwd=0.3)
          if (print.legend) {
            legend(max(gauss.variance), percent, 
                   xj=1, yj=1, ncol=2, 
                   legend=c("max", expression(l[2]), 
                     expression(l[1]), expression(r[1]), expression(r[2]), 
                     expression(r[2, 2]), expression(r[2, 3]),
                     expression(r[2, 4]), expression(r[2, 5]))[select.lp], 
                   col=col.leg, pch=15, cex=cex)
            if (type=="p")
              legend(max(gauss.variance),  y.pos.leg * percent, 
                     xj=1, yj=0.5, ncol=2, 
                     legend=c("max-min", "L2, no bin", "L1, no bin", 
                       expression(w[1]),  expression(w[2]), 
                       expression(w[3]), expression(w[4]), 
                       expression(w[5]), expression(w[6]), expression(w[7]), 
                       )[select.weight], 
                     pch=pch.leg[select.weight], cex=cex)
            else            
              legend(max(gauss.variance),  y.pos.leg * percent, 
                     xj=1, yj=0.5, ncol=2, 
                     legend=c("max-min", "L2, no bin", "L1, no bin", 
                       expression(w[1]),  expression(w[2]), 
                       expression(w[3]), expression(w[4]), 
                       expression(w[5]), expression(w[6]), expression(w[7]), 
                       )[select.weight], 
                     lty=lty.leg[select.weight], cex=cex)
          }
          s <- apply(X[[i]][[p]][without1, , drop=FALSE], 2, mean)
          r <- rank(s)
          if (PrintLevel>6) {
            if (i<=3) {## i.e., E1, V2, S3
              print(format(c(i, pvalue[p], 
                             s[c(1, .mpp.tests-1, .mpp.tests)], 
                             r[c(1, .mpp.tests-1, .mpp.tests)], Inf, 
                             mean(matrix(s[2:(.mpp.tests-2)], 
                                         nrow=.mpp.l.norms)) / number.files), 
                           dig=.mpp.digits), 
                    quote=FALSE)
              print(format(cbind(matrix(s[2:(.mpp.tests-2)], 
                                        nrow=.mpp.l.norms), 
                                 seq(Inf, Inf, l=.mpp.l.norms), 
                                 matrix(r[2:(.mpp.tests-2)], 
                                        nrow=.mpp.l.norms)), dig=.mpp.digits), 
                    quote=FALSE)
            } else { ##  MAX4
              print (cbind(.mpp.maxtests, s, r))
            }
          }
          Dev(FALSE, dev)
          if (is.numeric(dev) && !is.null(ps)) rl(filename)
        } ## p
      } ## i
    } ## print.data <- function(){...
    
    titleblock <- basicblock <- forblock <- forendblock <- "";
    for (i in 1:length(args)) {
      titleblock <-
        paste(titleblock, '"  ', names(args)[i], '=", args.', i, ', ', sep="")
      if (is.numeric(args[[i]][1])) {
        basicblock <- paste(basicblock, 'args.', i, '*100, ".", ', sep="")
      } else {
        basicblock <- paste(basicblock, 'args.', i, ', ".", ', sep="")
      }      
      forblock <-
        paste(forblock, 'for (args.', i, ' in args[[', i, ']]) {', sep="")
      forendblock <- paste(forendblock, "}", sep="")
    }
    titleblock <-
      paste('Titel <- paste(Title, ", ", ', titleblock, 'sep="");', sep="");
    if (PrintLevel>7) {print(titleblock);  print(basicblock);  print(forblock)}
    
    for (ind in npoints) { ##????
      BASICblock <- paste('basicname <- paste(name, ".", ', basicblock,
                          'ind, ".", sep="");', sep="")
      PSblock <- paste('psbasicname <- paste(ps, ".", ', basicblock,
                       'ind, sep="");', sep="")
      FORblock <- paste(forblock, titleblock, BASICblock, PSblock, 
                        "print.data();", forendblock, sep="");
      if (PrintLevel>7) print(FORblock)
      eval(parse(text=FORblock))
    }
    
    ntot <- 1
    for (i in 1:length(args)) ntot <- ntot * length(args[[i]])
    
    if (PrintLevel>7) print("TOTAL")
    ret <- array(dim=c(length(Numbers), length(pvalue), .mpp.tests))
    for (i in 1:length(Numbers)) {
      for (p in 1:length(pvalue)) {
        if (PrintLevel>7) print(paste(simu.fig1.ps.special[i], pvalue[p]))
        s <- 100 * apply(Xtotal[[i]][[p]][without1, , drop=FALSE], 2, sum) /
          (nrow(Xtotal[[i]][[p]][without1, , drop=FALSE]) * ntot)
        ret[i, p, ] <- s
        r <- rank(s)
        if (PrintLevel>7) {
          if (i<=3) {
            print(format(c(i, round(s[c(1, .mpp.tests - 1, .mpp.tests)]), 
                           r[c(1, .mpp.tests - 1, .mpp.tests)], Inf, 
                           sum(matrix(s[2:(.mpp.tests - 2)], 
                                      nrow=.mpp.l.norms)) ), 
                         dig=.mpp.digits), 
                  quote=FALSE)
            print(format(cbind(round(matrix(s[2:(.mpp.tests - 2)], 
                                            nrow=.mpp.l.norms)), 
                               seq(Inf, Inf, l=.mpp.l.norms), 
                               matrix(r[2:(.mpp.tests - 2)], 
                                      nrow=.mpp.l.norms)), 
                         dig=.mpp.digits), 
                  quote=FALSE)
          } else {
            print (cbind(.mpp.maxtests, s, r))
          }
        }
      }
  }
    return(ret)
  }

  
  
#### latex-table, e.g. table 1 in SRD (2004)
  x.table <- function(s, label, caption, reverse.ordering=FALSE, file, 
                      norm.select=c(1:5, 9), weight.select=1:.mpp.weights, 
                      parenthesis.space="4.12ex", 
                      add.col=NULL, start=TRUE, end=TRUE)
    {
      
      ## s : table of the rejection rates
      ## label : reference label of the latex table, "tab:" added
      ##         only used if start=TRUE
      ## caption : latex caption of table
      ## reverse.ordering: in parentheses behing the rejection rates, the
      ##          rank within the table is given. If FALSE the ranking is
      ##          given in the reverse ordering; this is of interest if
      ##          none of the simulations should be rejected
      ## file : filename of the latex file
      ## norm.select, weight.select : norms and weights that should be printed
      ## parenthesis.space : hbox for the ranking numbers
      ## add.col : if not NULL an additional column is inserted as first column
      ##           with add.col as entry in the first row 
      ## start, end : values that determine whether the table should be started
      ##              or termined
      col.names <- c("const", "$w_2$", "$w_3$", "$w_4$", "\\#", "$\\sqrt{\\#}$", 
                     "sd")[weight.select]
      col.names <- c("$w_1$", "$w_2$", "$w_3$", "$w_4$", "$w_5$", "$w_6$", 
                     "$w_7$")[weight.select]
      row.names <- c("$l_\\infty$", "$l_2$", "$l_1$", "$r_{1}$", "$r_{2}$", 
                     "$r_{{2}}$", "$r_{{{2}}}$", "$r_{{{{2}}}}$", 
                     "$r_{{{{{2}}}}}$")[norm.select]
      if (is.null(add.col))  addand <- ""
      else {
        addand <- "&"
      }
      stopifnot(!missing(file))
      write(file=file, append=!start, paste("% x.table for mpp.R <", add.col, 
              "> (start, end)=(", start, end, ")"))
      if (!missing(label) && start)
        label <- paste("\\label{tab:", label, "}", sep="")
      else label <- ""
      if (!missing(caption) && start) {
        write(file=file, append=TRUE,
              paste("\\caption{", caption, label, "}", sep=""))
      }
      if (start)
        write(file=file, append=TRUE, 
              paste("\\begin{tabular}{", 
                    paste(rep("r", length(weight.select)+1+!is.null(add.col)), 
                          collapse=""), 
                    "}"))
      s <- s[norm.select, weight.select]
      r <- round(matrix(rank( (-1 + 2 * reverse.ordering) * s),
                        ncol=.mpp.weights))
      s <- round(s)
      if (start) {
        write(file=file, append=TRUE, 
              paste(addand, "&", paste(col.names, collapse=" & "), "\\\\"))
        write(file=file, append=TRUE, "\\hline")
      }
      for (i in 1:length(norm.select)) {
        
        if (i==1) first.col <- add.col else first.col <- ""
        write(file=file, append=TRUE, 
              paste(first.col, addand, 
                    row.names[norm.select[i]], "&", 
                    paste(paste(s[i, ], " \\hbox to ", parenthesis.space, 
                                "{(\\hfill", r[i, ], ")}", sep=""), 
                          collapse=" & "), 
                    "\\\\"))
      }
      write(file=file, append=TRUE, "\\hline")     
      if (end) write(file=file, append=TRUE, "\\end{tabular}")
    }

  
######################################################################
######################################################################
  plotFig1 <- function(norm.select=simu.fig1.norm.select,
                       weight.select=simu.fig1.weight.select,
                       dev, coloured=!is.logical(dev),
                       individ = simu.fig1.individ,
                       absolute = TRUE,#if FALSE then alpha-value are based
                       ##             on the 100% Gauss model
                       estimate = TRUE, # estimate parameters by MLE or
                       ##             take the true ones?
                       ##             cf. function G in worksheet.add
                       simulate = TRUE, # simulate random coins? (if not
                       ##                take Gauss rf with same C)
                       ##             cf. function G in worksheet.add
                       MCrep = simu.MCrep
                       ) {
    pvalue <- 90
    bin <- c(-1, seq(0, 0.7, l=15))
    file.name <- simu.models$add$file
    Title <- simu.models$add$title
    mpp.model <- simu.models$add$model
    mar <- c(4.1, 4.1, 0.65, 0.4)  
    x <- sqrt(individ / 50)
    y <- sqrt(individ / 50)
    gauss.variance<- simu.gauss.variance
    
    
    psname <- paste(ps.path, file.name, sep="")
    xlim <- c(0, x)
    ylim <- c(0, y)
    if (estimate) {
      MCm <- list(model=simu.standard.estim, param=c(0, NA, 0, NA))
      sill <-  NA
    } else {
      MCm <- model
      sill <- 1
    }
     
    for (r in 1:length(simu.rlist)) {
      if (PrintLevel>1) cat("plot.fig1", r, "\n")
      for (leg in c(FALSE, TRUE)) {
        if (PrintLevel>3) print(leg)
        par(mar=mar)
        plot.fig1(name=paste(simu.path, file.name, sep=""), 
                  npoints=individ, 
                  gauss.variance=simu.gauss.variance, 
                  Title=Title,
                  ## pvalue,
                  ps = paste(psname, leg, sep="."),
                  absolute=absolute, 
                  ## PrintLevel
                  select.lp=norm.select, 
                  select.weight=weight.select, 
                  ##select.extra, return.rate
                  dev = dev,
                  ##print.title
                  height=4.5, width=5,                     
                  coloured=coloured,
                  print.legend=leg, 
                  type="b", cex=fig1.cex,
                  y.pos.leg =0.55,
                  pch.basic=c(1, 0, 4, 0, 0),
                  lty.basic=c(0, 0, 2, 0, 1, 0, 3), 
                  lty.abline=4,
                  ## col.abline
                  
                  estimate=estimate, 
                  simulate=simulate, 
                  MCmodel=simu.standard.estim,
                  gauss.model=simu.standard.cov,
                  radius=simu.rlist[[r]]$radius, 
                  scale=simu.rlist[[r]]$scale
                  )
        par(mar=mar)
      }
    } # r in simu.rlist

    simu <- list()
    if (simu.rlist[[r]]$simu) {
      for (r in 1:length(simu.rlist)) {
        if (PrintLevel>8) cat("simulate.rfmtest", r, "\n")
        model <- function(variance) # 
          list(list(model = mpp.model, 
                    p = c(fct=1, # coin function
                      radius=simu.rlist[[r]]$radius, 
                      height=sqrt((1-variance) * x *
                        y / (individ*pi*simu.rlist[[r]]$radius^2)))), 
                      "+", 
               list(model=simu.standard.cov, var=variance,
                    scale=simu.rlist[[r]]$scale)
               )
        simu[[r]] <-
          simulate.rfmtest(coord=NULL,  npoints=individ, 
                           coordmodel="uniform",  window=c(xlim, ylim), 
                           repetitions=1, edgecorrection=simu.rlist[[r]]$radius,
                           coordrep=1, model = model, 
                           gauss.variance = simu.rlist[[r]]$gauss.variance,
                           normalize=TRUE, MCrepetition=MCrep, MCmodel = MCm, 
                           sill= sill, bin=bin, Ebin=seq(0, 1, 0.01), 
                           saving=FALSE, path=simu.path, return.data=TRUE)
      }

      for (r in 1:length(simu.rlist)) {
        filename <- paste(psname, "simu.circles", simu.rlist[[r]]$radius * 100,
                          simu.rlist[[r]]$scale * 100, 
                          simu.rlist[[r]]$gauss.variance, sep=".")
        Dev(TRUE, dev, ps=filename)
        plotWithCircles(cbind(simu[[r]]$coord, 
                              0.3+simu[[r]]$data - min(simu[[r]]$data)), 
                        xlim=xlim, ylim=ylim, factor=0.01) # factor
        Dev(FALSE, dev)
        if (is.numeric(dev)) rl(filename)
      }
    }
  }


  ######################################################################
  data.basics <- function(data = NULL, coord=NULL, 
                          radius=NULL, scale=0, sill=0, individ=50, 
                          simuxlim =c(0, 1), simuylim = c(0, 1), 
                          bin=data.bin, 
                          MCmodel = list(model="exponential", 
                            param = c(0, NaN, 0, NaN)), 
                          MCrep = simu.MCrep, 
                          lower.kappa = NULL, upper.kappa = NULL, 
                          normalize=TRUE, showEfct = FALSE,
                          name=paste("dummy", radius*100, scale*100, individ, 
                            sill*100, nr, sep="."), 
                          ## plotting with radii
                          height=5, 
                          pfactor=diff(simuxlim)*diff(simuylim)
                          ) {
    ## dev = 2,  ## global
    par(parlist)
    if (PrintLevel>3) cat(name, "\n")
    ps.path.name <- paste(ps.path, name, sep="")
    filename <- paste(data.save.path, 
                      name, ".", paste(c(MCmodel, recursive=TRUE), collapse="."),
                      ".data", sep="")
    if (file.exists(filename)) {
      if (PrintLevel>0) cat("using values stored in", filename, "\n")
       load(filename)
    } else {
      stopifnot(!is.null(data))
      x <- list()
      if (is.null(coord)) {
        x$coord <- data[, c(1, 2)]                       
        x$data <- data[, c(-1, -2), drop=F]
      } else {
        x$coord <- coord
        x$data <- as.matrix(data)
      }
    }

    x$NormedData <- apply(x$data[, 1, drop=F], 2, norm)## only a good approximation
    ##                                          if data are uniformly scattered
  
    ## plot the data
    width <- height * diff(range(x$coord[, 1])) / diff(range(x$coord[, 2]))
    Dev(TRUE, dev, ps=ps.path.name, height=height, width=width)
    cat("file name ", ps.path.name, "", dev, width, "\n")
    if (FALSE) {
      r <- x$NormedData - min(x$NormedData) + 0.5
      r <- r * pfactor / (max(r) * nrow(x$coord));
    } else {
      r <- x$data * pfactor
    }
    plotWithCircles(cbind(x$coord, r), 
                    xlim=range(x$coord[, 1]), ylim=range(x$coord[, 2])
                    )
    Dev(FALSE, dev)
    if (is.numeric(dev)) rl(filename)

     
    ## rfm.test
    if (!file.exists(filename)) {
      if (PrintLevel>6) cat("rfm.test")
      y <- rfm.test(# coord
                    data=x, normalize=normalize, MCrepetition=MCrep, 
                    MCmodel=MCmodel,
                    # method,
                    bin=bin, Ebin=seq(0, 1, 0.01),
                    # use.naturalscaling, MCregister,
                    n.hypo=data.n.hypo,
                    # pvalue, tests, tests.lp, tests.weight, Barnard, Print
                    optim.control = optim.control,
                    lower=lower.kappa, upper=upper.kappa
                    )
      if (!debug) y$null.hypo <- NULL
      save(file=filename, x, y)
      if (PrintLevel>9) str(y)
    }
    
    ## show characteristics
    est <- y$est[[1]]
    filename <- paste(filename, "2", sep="")
    if (file.exists(filename)) load(filename)
    if (!file.exists(filename) || showEfct) {
      if (PrintLevel>6) cat(" mpp.characteristics")
      z <- mpp.characteristics(bin=bin, rep=1, p=0.8, 
                               name=ps.path.name,
                               show=showEfct,
                               normalize=normalize, model=est, 
                               c1=x$coord, d1=x$data,
                               readline=rl,
                               dev=dev)
    }
    if (!file.exists(filename)) {
      if (PrintLevel>6) cat("\nVariogram analysis\n")
      zz <- GaussRF(x$coord, grid=FALSE, model=est, n=data.conf.repet)
      b <- matrix(ncol=data.conf.repet, nrow=length(bin)-1)
      for (i in 1:data.conf.repet)
        b[,i] <- EmpiricalVariogram(x=x$coord, data=norm(zz[, i]),
                                    grid=FALSE, bin=bin)$emp
      save(file=filename, b, z)
    }
    if (PrintLevel>6) cat("\n\n")

    v <- Variogram(z$midbin, model=est)
    v[1] <- NA ## do not plot for 0
    lower <- apply(b, 1, quantile, data.confid.p / 2)
    lower[1] <- NA
    upper <- apply(b, 1, quantile, 1 - data.confid.p / 2)
    upper[1] <- NA
    Dev(TRUE, dev, ps=paste(ps.path,name,".confidence.gamma",sep=""),
        height=3)
    matplot(z$midbin,
            cbind(lower , v, upper, z$GAM),
            col=c("black","black","black","black"), ylim=c(0.5, 1.1),
            pch=c(NA, NA, NA, 16), lty=c(2, 1, 2, 0),
            type=c("l", "l", "l", "p"), xlab="distance", ylab=expression(gamma))
    Dev(FALSE, dev)
    
    return(list(data=x, analyse=y, plots=z, b=b))
  } # data.basics
  
  plotdata <- function(data) {
    par(parlist)
    ALL <- list()
    ## name="steiger"
    ## name="coul"
    for (i in 1:length(data)) {
      z <- data.basics(data=get(data[[i]]$set)$diameter, 
                       coord=get(data[[i]]$set)$coord, 
                       showEfct=get(paste(data[[i]]$set, ".showEfct", sep="")), 
                       MCmodel=data.MCmodel, 
                       lower.kappa = 0.01, upper.kappa=10, 
                       normalize=get(paste(data[[i]]$set, ".showEfct", sep="")), 
                       name=data[[i]]$set, 
                       pfactor=data[[i]]$pfactor ## 2: diameter!
                       )
      ALL <- append(ALL, 
                    list(list(d=z, lty=i, lwd=2, col="black",
                              pch=c(2, 16, 3)[i], legend=data[[i]]$name
                              )))
    }
    
    leg <- NULL
    for (f in 1:length(data.fcts)) {
      if (PrintLevel>4) cat("\n", data.fcts[f], "...")
      filename <- paste(ps.path, "data.", data.fctnames[[f]], sep="")
      Dev(TRUE, dev, ps=filename, height=4, width=3.5)
      plot(Inf, Inf, xlim=data.xlim, ylim=data.ylims[, f], xlab=data.xlab, 
           ylab=data.ylabs[f])
      for (i in 1:length(ALL)) {
        eval(parse(text=paste("v <- ALL[[i]]$d$plots$", data.fcts[f], sep="")))
        if (PrintLevel>7) {
          print(as.vector(v))
          print(ALL[[i]]$d$plots$midbin)
        }
        points(ALL[[i]]$d$plots$midbin, v, 
               pch=ALL[[i]]$pch, col=ALL[[i]]$col, 
               lty=ALL[[i]]$lty, lwd=ALL[[i]]$lwd)
        if (data.fcts[f]=="GAM") {
          leg <- c(leg, ALL[[i]]$legend)
          x <- ALL[[i]]$d$analyse$bin
          x <- x[x>0]
          lines(x, Variogram(x, model=ALL[[i]]$d$analyse$est[[1]]), 
                lty=ALL[[i]]$lty, lwd=ALL[[i]]$lwd, col="black") # col="black"
        }
      }
      if (data.fcts[f]=="GAM")
        legend(x=max(data.xlim), y=0, xj=1, yj=0, legend=leg, 
               lty=sapply(ALL, function(x) x$lty), 
               cex=1.4)
      if (is.numeric(dev)) rl(filename)
      Dev(FALSE, dev)
    }
  } # plotdata


  ######################################################################
  theoretical.examples <- function() {
    par(parlist)    
    coin.E <- function(x, lambda=1, R=1, d=2) {
      bd <- pi^(d/2) / gamma(1 + d/2)
      res <- rep(1 + lambda * bd * R^d, length(x))
      sel <- (x!=0) & (abs(x)<=R)
      res[sel] <- res[sel] + 1.0
      return(res)
    }
    coin.V <- function(x, lambda=1, R=1, d=2) {
      bd <- pi^(d/2) / gamma(1 + d/2)
      lambda * bd * R^d
    }
    coin.G <- function(x, lambda=1, R=1, d=2) {
      stopifnot(d==2)
      bd <- pi^(d/2) / gamma(1 + d/2)
      x <- abs(x / (2 * R))
      res <- rep(1, length(x))
      res[x<1] <- 1 - 2 / pi * (acos(x[x<1]) - x[x<1] * sqrt(1- x[x<1]^2))
      return(res * (lambda * R^d * pi))
    }

    lambda <- 1
    R <- 1 
    x <- seq(0, 2.8, len=theo.len)
    Dev(TRUE, dev, ps=paste(ps.path, "coin.EVG", sep=""), height=5)
    matplot(x, cbind(coin.V(x=x, lambda=lambda, R=R), 
                     coin.G(x=x, lambda=lambda, R=R)            
                     ), 
            lty=c(theo.ltyV, theo.ltyG), type="l", lwd=theo.lwd, col=1, 
            pch=theo.pch, xlab="", ylab="", ylim=c(0, 5.1), cex.axis=2.3,
            cex=1.5)
    points(0, coin.E(0), pch=theo.pchE)
    points(0, coin.V(0), pch=theo.pchV)
    points(0, coin.G(0), pch=theo.pchG)
    x1 <- seq(0.000001, 0.9999, len=50)
    lines(x1, coin.E(x=x1, lambda=lambda, R=R), lty=theo.ltyE)
    x1 <- seq(1.000001, 2.8, len=50)
    lines(x1, coin.E(x=x1, lambda=lambda, R=R), lty=theo.ltyE)
    Dev(FALSE, dev)
    
    var.E <- function(x, lambda=1, R=1, d=2) rep(0, length(x))
    var.V <- function(x, lambda=1, R=1, d=2) {
      bd <- pi^(d/2) / gamma(1 + d/2)
      res <- rep(1 + lambda * bd *R^d * (3 + lambda * bd * R^d), length(x))
      sel <- (x!=0) & (abs(x)<=R)
      res[sel] <- res[sel] + 3 + 2 * lambda * bd * R^d
      return(res)
    }
    var.G <- function(x, lambda=1, R=1, d=2) as.integer(x!=0) *
      var.V(0, lambda=lambda, R=R, d=d)
    
    lambda <- 1
    R <- 1 
    x <- seq(0.00001, 2.1, len=theo.len)
    Dev(TRUE, dev, ps=paste(ps.path, "var.EVG", sep=""), height=5)
    matplot(x, cbind(var.E(x=x, lambda=lambda, R=R), 
                     var.G(x=x, lambda=lambda, R=R)           
                     ), 
            lty=c(theo.ltyE, theo.ltyG), 
            type="l", lwd=theo.lwd, col=1, pch=theo.pch, xlab="", ylab="", 
            cex.axis=2.3, cex=1.5, ylim=c(0, 30)
            )
    points(0, var.E(0), pch=theo.pchE)
    points(0, var.V(0), pch=theo.pchV)
    points(0, var.G(0), pch=theo.pchG)
    x1 <- seq(0.000001, 0.9999, len=50)
    lines(x1, var.V(x=x1, lambda=lambda, R=R), lty=theo.ltyV)
    x1 <- seq(1.000001, 2.8, len=50)
    lines(x1, var.V(x=x1, lambda=lambda, R=R), lty=theo.ltyV)
    Dev(FALSE, dev)
    
    nn.E <- function(x, lambda) {
      res <-  rep( 1 / (2 * sqrt(lambda)), length(x))
      for (i in 1:length(x)) if (x[i]!=0) {
        res[i] <-
          integrate(f=function(s) exp(-s)/sqrt(s), 
                    low=1E-100, up=lambda * pi * x[i]^2)$value /
                      (2 * sqrt(pi * lambda))
      }
      return(res)
    }
    nn.V <- function(x, lambda) (1 - exp(-lambda * pi * x^2) + (x==0)) /
      (lambda * pi) - nn.E(x, lambda)^2
    
    nn.G <- function(x, lambda) {
      U <- function(r, s, t, d=2) {
        stopifnot(all(c(r, s, t)>=0))
        bd <- pi^(d/2) / gamma(1 + d/2)
        l <-  max(c(length(r), length(s), length(t)))
        r <- rep(r, length=l)
        s <- cbind(rep(s, length=l), t)
        res <- rep(0, nrow(s))
        r0 <- r==0
        if (any(r0)) {
          res[r0] <- bd * apply(s[r0, , drop=FALSE], 1, max)^d
        }
        if (!all(r0)) {
          for (i in 1:2) {
            s1 <- s[, i]
            s2 <- s[, 3 - i]
            x <- s1
            sel1 <- (s1 <= s2) & (s1 <= s2 - r)
            x[sel1] <- -s1[sel1]
            sel <- !sel1 & ( (s1 + s2)/r - 1 > 1E-12) & ((s1 < s2) | (s1 < s2+r))
            x[sel] <- (s1[sel]^2 + r[sel]^2 - s2[sel]^2) / (2 * r[sel])
            stopifnot(d==2)
            
            if (x^2 > s1^2) {
              print(x); print(s1); print(x-s1); print(s2); print(r);
              print(sel); print(sel1); print(s1<=s2); print(s1 + s2 - r);
              stop(" @@ ")
            }
            S <- x * sqrt(s1^2 -x^2) + s1^2 * asin(x/s1) + s1^2 * pi / 2
            res[!r0] <- res[!r0] + S[!r0]
          }
        }
        return(res)
      }
      gam <- (1 - exp(-lambda * pi * x^2) ) / (lambda * pi)
      for (i in 1:length(x)) {
        a <- adapt(2, lower=c(0, 0), upper=c(x[i], x[i]), min=100, max=1000, 
                   fun=function(st) exp(-lambda * U(x[i], st[1], st[2])), 
                   eps=0.01)
        gam[i] <- gam[i] - a$value
      }
      return(gam)
    }
    
    x <- seq(0.0001, 1.8, len=theo.len)
    lambda <- 1
    R <- 1
    if (!exists("nnE")) nnE <- nn.E(x=x, lambda=lambda)
    if (!exists("nnV")) nnV <- nn.V(x=x, lambda=lambda)
    if (!exists("nnG")) nnG <- nn.G(x=x, lambda=lambda)
    Dev(TRUE, dev, ps=paste(ps.path, "nn.EVG", sep=""), height=5)
    matplot(x[-1], cbind(nnE[-1], nnV[-1], nnG[-1]), 
            lty=c(theo.ltyE, theo.ltyV, theo.ltyG), type="l", lwd=theo.lwd,
            col=1, pch=theo.pch, cex.axis=2.3, cex=1.5, xlab="", ylab="")
    points(0, nn.E(0, lambda=lambda), pch=theo.pchE, cex=1.5)
    points(0, nn.V(0, lambda=lambda), pch=theo.pchV, cex=1.5)
    points(0, nn.G(0, lambda=lambda), pch=theo.pchG, cex=1.5)
    Dev(FALSE, dev)
  }

  #########################################

  worksheet.add <- function(i=50,           # intensity
                            x=sqrt(i / 50), # length of window
                            y=sqrt(i / 50), # width of window
                            w=500,          # number of repetitions
                            model = simu.models$add,
                            Barnard = FALSE, # application of test by
                            ##  Barnard & Besag/Diggle, instead of that by SRD
                            R.vario = simu.R.vario,
                            radii = simu.radii
                            ) {
    mpp.model <- model$model
    basic.name <- model$file
    G <- function(R=0.1, w=1000, i=50, x=1, y=1, simulate=TRUE, estimate=TRUE, 
                  radius=0.1, MCrep=simu.MCrep, MCmodel=simu.standard.MC,
                  gauss.model=simu.standard.cov,
                  method=if (x * y * i > 500 || gauss.model=="wave" ||
                    MCmodel=="wave")
                  NULL else "direct"){      
      if (PrintLevel>1)
        cat("G ", gauss.model, "(", MCmodel, ")", " e=", estimate,
            " s=", simulate, " radius=", radius, " R=", R,
            " intens.=", i, " len=", x, " wid=", y, "\n",
            sep="")
      ## R: scale of variogram
      ## w: number of replicat.s for each parameter set (should ideally be 5000)
      ## i: number of individuals (points/trees) in each simulation
      scale <- R
      rep <- as.integer(w)
      individ <- as.integer(i)
      xlim <- c(0, x)
      ylim <- c(0, y)
      names <- GetModelNames()
      MCmodel <- names[pmatch(MCmodel, names)]
      gauss.model <- names[pmatch(gauss.model, names)]
      
      name <- simu.addname(basic.name, estimate, simulate, MCmodel, gauss.model,
                           radius, scale, individ)
      
      if (simulate) {
        lbRd <- individ * pi * radius^2 / (diff(xlim) * diff(ylim))
        intrinsic.var <- model$intrinsic.var(lbRd)
        model <- function(variance)
          list(list(model=mpp.model, 
                    p=c(fct=1, # coin function
                      radius=radius, 
                      height=sqrt((1-variance)  / intrinsic.var))), 
               "+", list(model=gauss.model, var=variance, scale=scale))
      } else {
        ## the local variables of scale and radius are used here!
        model <- function(variance)
          list(list(model="circ", var=1.0-variance, scale=2 * radius), 
               "+", list(model=gauss.model, var=variance, scale=scale))
      }
      
      if (estimate) {
        MCm <- list(model=MCmodel, param=c(0, NA, 0, NA))
        sill <-  NA
      } else {
        ## scale is really the local scale here !
        MCm <- function(variance)
          list(model=gauss.model, param=c(0, variance, 0, scale))
        sill <- NA
      }

      simulate.rfmtest(npoints=individ, 
                       window=c(xlim, ylim), edgecorrection=radius, 
                       coordrep=rep, model = model,
                       normalize=TRUE,
                       method=method,
                       MCrepetition=MCrep, MCmodel = MCm, sill= sill,
                       saving=name, path = simu.path,
                       Barnard = Barnard
                       )
    }
    
   if (!Barnard) {
      for (R in c(0.0001, R.vario))
        for (radius in radii) {
          G(e=TRUE, s=TRUE, radius=radius, R=R, w=w, i=i, x=x, y=y)
          GC <- gc()
          if (PrintLevel>6) print(GC)
        }
    }
    
    if (mpp.model==simu.models$add$model) {
      ##  Gaussian rf with variogram, estimated by wrong variogram model;
      ##  estimated parameters;   #simulate=FALSE
      for (R in R.vario)
        for (radius in radii) {
          G(e=TRUE, s=FALSE, radius=radius, R=R, w=w, i=i, x=x, y=y)
          if (PrintLevel>6) print(gc())
        }
    }
    
    if (mpp.model==simu.models$add$model  && !Barnard) {
      ## estimation of Gaussian rf with same variogram as additive boolean
      ## plus Gaussian random field. -- true parameters ##
      for (R in R.vario)
        for (radius in radii) {
          G(e=FALSE, s=FALSE, radius=radius, R=R, w=w, i=i, x=x, y=y)
          GC <- gc()
          if (PrintLevel>6) print(GC)
        }
      
      ## do not use complicated variograms. parameters cannot be transfered
      ## by est.gauss.add.work
      for (m in  list(c(gauss.m="cubic", MCm="expon"),
                      c(gauss.m="gneit", MCm="expon"),
                      c(gauss.m="gneit", MCm="wave"),
                      c(gauss.m="wave",  MCm="gneit"),
                      c(gauss.m="expon", MCm="gneit")))
        for (R in R.vario) {
          G(e=TRUE, s=FALSE, gauss.m=m[1], MCm=m[2], R=R, w=w, i=i, x=x,y=y)
          GC <- gc()
          if (PrintLevel>6) print(GC)
        }
    } # mpp.model=="random coin"
  } # worksheet.add


  worksheet.nn <- function(i=50,           # intensity
                           x=sqrt(i / 50), # length of window
                           y=sqrt(i / 50), # width of window

                           w=500,          # number of repetitions
                           Barnard = FALSE # application of test by
                           ##  Barnard & Besag/Diggle, instead of that by SRD
                           ) {
    
    Gnn <- function(R=0.1, w=1000, i=50, x=1, y=1, 
                    radius=0.1, gaussvariance=simu.gauss.variance,
                    MCrep=simu.MCrep, MCmodel=simu.standard.MC,
                    gauss.model=simu.standard.cov){
      if (PrintLevel>1) cat("Gnn R=", R, " repet=", w, " intensity=", i,
                            " length=", x, " width=", y, "\n", sep="")
      scale <- R
      rep <- as.integer(w)
      individ <- as.integer(i)
      xlim <- c(0, x)
      ylim <- c(0, y)
      name <- simu.nnname(scale, x, y, individ)
      
      model <- function(variance)
        list(list(model=gauss.model, var=variance, scale=scale), 
             "+", 
             list(model="nearest neighbour", 
              p=c(factor=sqrt((1 - variance) *
                    individ / (diff(xlim) * diff(ylim)) / (1/pi - 0.25))))
             )
      simulate.rfmtest(npoints=individ, window=c(xlim, ylim), 
                       edgecorrection=radius, coordrep=rep, 
                       model = model,
                       normalize=TRUE, 
                       MCrepetition=MCrep, 
                       MCmodel = list(model=MCmodel, param=c(0,NA,0,NA)), 
                       sill= NA,
                       saving=name, path = simu.path,
                       Barnard = Barnard
                       )
    }
    
    for (p in simu.R.nn.point.setting)
      for (R in  simu.R.nn.vario) {
        Gnn(R=R, w=w, i=p[1], x=p[2], y=p[3])
        GC <- gc()
        if (PrintLevel>6) print(GC)
      }
  } # worksheet.nn  

  latex.tab.123 <- function(force=FALSE,  # dev,
                            norm.select=1:3, #norm.select <- c(1:5, 9)
                            weight.select=1:.mpp.weights,
                            EVS = c(1, 2) #1:E, 2:V, 3:S
                            ) {
    lname <-  c("1", "2", "3", "4", "5", "6", "10", "11")
    lname <- lname[c(1:6)] ## print the simuluations that match
    ##                        the value of `label' in the list below
    Pvalues <- c(90, 95)
    pvalues <- c(2) #1: 90, 2: 95
    par(parlist)
    
    X <-
      list(list(file.name = simu.models$add$file, 
                label = "1", 
                Title = "Random coin model of Example 1", 
                tex.label = "add", 
                caption = paste("The average is over all combinations of radii",
                  "$R=", paste(simu.radii, collapse=","),
                  "$, $\\alpha=0, 0.1, \\ldots, 0.9$ and scales $s=",
                  paste(simu.R.vario, collapse=","),
                  "$ of the covariance function."), 
                null.hypo.text = 2, 
                variable =
                c(list(list(individ=100, EVS=1, add.col="$E$, $n=100$",
                          null.hypo=0)),
                  if (200 %in% simu.individ.list)
                  list(list(individ=200, add.col="\\hphantom{$E$, } $n=200$")), 
                  if (50 %in% simu.individ.list)
                  list(list(individ=50, add.col="\\hphantom{$E$, } $n=50$")), 
                  list(list(individ=100, EVS=2, add.col="$V$, $n=100$",
                            null.hypo=1)),
                  if (200 %in% simu.individ.list)
                  list(list(individ=200, add.col="\\hphantom{$V$, } $n=200$")), 
                  if (50 %in% simu.individ.list)
                  list(list(individ=50, add.col="\\hphantom{$V$, } $n=50$")), 
                  ),
                ),

           list(file.name = simu.models$nn$file, 
                label = "2", 
                Title = paste("Nearest neighbour model of Example 3. As in",
                  "Table \\ref{tab:add.1} the average of the rejection",
                  "rates are given, with ranks in parentheses. "),
                tex.label = "nn", 
                no.further.text=TRUE, 
                caption = paste("Simulations were of $n=100$ points on squares",
                  "with area $2$, $4$, and $8$; the scale parameter of the",
                  "random field took the values $s=",
                  paste(simu.R.nn.vario, collapse=","), "$."),
                individ = 100,
                scale = simu.R.nn.vario, 
                area = c(2, 4, 8)/100, 
                null.hypo= 0, 
                variable = list(list(EVS=1, add.col="$E$"), 
                  list(EVS=2, add.col="$V$")), 
                ), 
           
           list(file.name = simu.models$var$file, 
                label = "3", 
                Title = paste("Model given by Example 2. See Table",
                  "\\ref{tab:add.1} for the parameter choices and further",
                  "comments. "), 
                tex.label = "var", 
                no.further.text = TRUE, 
                caption="",#"See Table \ref{tab:add.1} for the parameter.", 
                null.hypo.text = 1, 
                variable =
                list(list(individ=100, EVS=1, add.col="$E$, $n=100$",
                          null.hypo=1), 
                     list(individ=100, EVS=2, add.col="$V$, $n=100$",
                          null.hypo=0), 
                     ), 
                ), 
           
           list(file.name = simu.models$add$file, 
                label = "4", 
                Title = paste("Gaussian random field model with the same",
                  "covariance structure as the random coin model"), 
                tex.label = "add", 
                caption = paste("In the MC test the exponential model is",
                  "fitted. (a)  $R=", paste(simu.radii, collapse=","),
                  "$ and exponential covariance structure of the random field",
                  "with $s=", paste(simu.R.vario[-1], collapse=","),
                  "$. (b)  $R=0.1$ and cubic or Gneiting's covariance structure",
                  "with  $s=", paste(simu.R.vario[-1], collapse=","), "$."), 
                absolute = TRUE, 
                simulate = FALSE, 
                scale = simu.R.vario[-1], 
                null.hypo.text = c(1, 2), 
                null.hypo= 1, 
                
                variable =
                c(if (50 %in% simu.individ.list)
                  list(list(EVS = 1, 
                            individ=50, radius=simu.radii, 
                            simu.cov = simu.standard.cov, 
                            add.col="$E$ (a) $n=50$")), 
                  list(list(EVS = 1,
                            individ=100, radius=simu.radii, 
                            simu.cov = simu.standard.cov, 
                            add.col="\\hphantom{(a) }$n=100$")), 
                  if (50 %in% simu.individ.list)
                  list(list(individ=50, radius = c(0.10), 
                       simu.cov = c("cubic", "gneiting"), 
                       add.col="$E$ (b) $n=50$")), 
                  list(list(individ=100, add.col="\\hphantom{(b) }$n=100$")), 
                  if (50 %in% simu.individ.list)
                  list(list(EVS = 2, 
                       individ=50, radius=simu.radii,
                       simu.cov = simu.standard.cov, 
                       add.col="$V$ (a) $n=50$")), 
                  list(list(individ=100, add.col="\\hphantom{(a) }$n=100$")), 
                  if (50 %in% simu.individ.list)
                  list(list(individ=50, radius = c(0.10), 
                       simu.cov = c("cubic", "gneiting"), 
                       add.col="$V$ (b) $n=50$")), 
                  list(list(individ=100, add.col="\\hphantom{(b) }$n=100$"))
                  ),           
                ), 
                 
           list(file.name = simu.models$add$file, 
                label = "5", 
                Title = "Random coin model", 
                tex.label = "add", 
                caption = paste("The average is over all combinations of $",
                  "R=0.05, 0.1$** and $s=", paste(simu.R.vario, collapse=","),
                  "$ and exponential covariance structure of the random field.",
                  "In th MC test an exponential model is fitted."), 
                radius = c(0.05, 0.10), 
                null.hypo=0, 
                EVS = 1, 
                variable =
                c(if (50 %in% simu.individ.list)
                  list(list(individ = 50, add.col = "n=50")), 
                  list(list(individ = 100, add.col = "n=100")), )
                ), 
           
           list(file.name = simu.models$add$file, 
                label = "6", 
                Title = "Random coin model", 
                tex.label = "add", 
                caption = paste("The average is over all combinations of **$",
                  "R=0.05, 0.1$** and $s=", paste(simu.R.vario, collapse=","),
                  "$ and exponential covariance structure of the random field.",
                  "In th MC test an exponential model is fitted."), 
                radius =  c(0.05, 0.1),
                null.hypo=0, 
                EVS = 1, 
                variable =
                c(if (50 %in% simu.individ.list)
                  list(list(individ = 50, add.col = "$n=50$")), 
                  list(list(individ = 100, add.col = "$n=100$")), )
                ), 
           list(file.name = simu.models$add$file,  
                label = "10", 
                Title = paste("Gaussian random field model with the same",
                  "covariance structure as the random coin model"), 
                tex.label = "add", 
                caption = paste("In the MC test the exponential model is",
                  "fitted. (a)  $R=", paste(simu.radii, collapse=","),
                  "$ and GNEITING covariance structure of the random field with",
                  "$s=", paste(simu.R.vario[-1], collapse=","),
                  "$. (b)  $R=0.1$ and cubic or Gneiting's covariance structure",
                  "with  $s=", paste(simu.R.vario[-1], collapse=","), "$."), 
                absolute = TRUE, 
                simulate = FALSE, 
                scale = simu.R.vario[-1], 
                null.hypo= 1, 
                radius = c(0.10), 
                
                variable =
                c(if (50 %in% simu.individ.list)
                  list(list(EVS = 1, 
                            individ=50, 
                            estim.cov = c("gneiting"), 
                            simu.cov = simu.standard.cov, 
                            add.col="E (e/g) n=50")), 
                  list(list(EVS = 1,
                            individ=100,
                            estim.cov = c("gneiting"), 
                            simu.cov = simu.standard.cov,
                            add.col="\\hphantom{(a) }n=100")), 
                  if (50 %in% simu.individ.list)
                  list(list(individ=50, 
                            simu.cov = c("wave"), 
                            add.col="E (w/g) n=50")), 
                  list(list(individ=100, add.col="\\hphantom{(b) }n=100")), 
                  if (50 %in% simu.individ.list)  
                  list(list(individ=50, 
                            simu.cov = simu.standard.cov, 
                            estim.cov = simu.standard.estim, 
                            add.col="E (e/e) n=50")), 
                  list(list(individ=100, add.col="\\hphantom{(b) }n=100")),   
                  if (50 %in% simu.individ.list)
                  list(list(EVS = 1, 
                            individ=50, 
                            estim.cov = c("wave"), 
                            simu.cov = c("gneiting"), 
                            add.col="E (g/w) n=50")), 
                  list(list(individ=100, add.col="\\hphantom{(a) }n=100")), 
                  )
                ), 
                
           list(file.name = simu.models$add$file,  
                label = "11", 
                Title = paste("Gaussian random field model with the same",
                  "covariance structure as the random coin model"),
               tex.label = "add", 
                caption = paste("In the MC test the", simu.standard.MC,
                  "model is fitted. (a)  $R=", paste(simu.radii, collapse=","),
                  "$ and GNEITING covariance structure of the random field with",
                  "$s=", paste(simu.R.vario[-1], collapse=","), "$. (b) $R=0.1$",
                  "and cubic or Gneiting's covariance structure with  $s=",
                  paste(simu.R.vario[-1], collapse=","), "$."),
                absolute = TRUE, 
                simulate = FALSE, 
                scale = simu.R.vario[-1], 
                null.hypo= 1, 
                radius = c(0.10), 
                
                variable =
                c(if (50 %in% simu.individ.list)
                  list(list(EVS = 2, 
                          individ=50, 
                          estim.cov = c("gneiting"), 
                          simu.cov = simu.standard.cov, 
                          add.col="V (e/g) n=50")), 
                  list(individ=100, add.col="\\hphantom{(a) }n=100"), 
                  if (50 %in% simu.individ.list)
                  list(list(individ=50, 
                       simu.cov = c("wave"), 
                       add.col="V (w/g) n=50")), 
                  list(list(individ=100, add.col="\\hphantom{(b) }n=100")), 
                  if (50 %in% simu.individ.list)  
                  list(list(individ=50, 
                            simu.cov = simu.standard.cov, 
                            estim.cov = simu.standard.estim, 
                            add.col="V (e/e) n=50")), 
                  list(list(individ=100, add.col="\\hphantom{(b) }n=100")), 
                  if (50 %in% simu.individ.list)
                  list(list(EVS = 1, 
                        individ=50, 
                          estim.cov = c("wave"), 
                          simu.cov = c("gneiting"), 
                          add.col="V (g/w) n=50")), 
                  list(list(individ=100, add.col="\\hphantom{(c) }n=100")), 
                  ), 
                )
           )
    
    for (p in pvalues) {
      for (ln in 1:length(lname)) {

        for (i in 1:length(X)) {        
          if (!is.na(pmatch(lname[ln], X[[i]]$label))) break
        }
        stopifnot(!is.na(pmatch(lname[ln], X[[i]]$label)))
        Y <- X[[i]]
        if (PrintLevel>5) {
          cat("\ni=", i, " ")
          if (PrintLevel>7) str(Y)
        }
        if (is.null(Y$variable)) Y$variable <- list(list(dummy=0))
        EVSgiven <- !is.null(Y$EVS)
        
        m <- list()
        file.name <-
          paste(simu.path, Y$file.name, ".", Y$label, ".analyse.summary.dat",
                sep="")
        if (file.exists(file.name)) load(file.name)
        tex.file.name <-
          paste(tex.path, Y$file.name, ".", Y$label, ".tex", sep="")  #
        first.null.hypo <- TRUE
        for (v in 1:length(Y$variable)) {
          cat("\n", p, "", ln, "", i, "", v)
          variable <- Y$variable[[v]]
          var.names <- names(variable)
          for (w in 1:length(variable)) {
            if (PrintLevel>2) cat(w, var.names[w], names(Y), v, "\n")
            stopifnot( (sum(var.names[w]==names(Y))==(v!=1)) )
            txt <- paste("Y$", var.names[w], "<-variable$", var.names[w], sep="")
            if (PrintLevel>5) cat(txt, "\n")
            eval(parse(text = txt))
          }
          if (!is.null(Y$gauss.variance)) gauss.variance <- Y$gauss.variance
          else gauss.variance <- simu.gauss.variance
          
          if (!file.exists(file.name)) {
            ## plot.fig1 is used only for the returned data
            m[[v]] <-
              plot.fig1(name=paste(simu.path, Y$file.name, sep=""), 
                        npoints=Y$individ, 
                        gauss.variance=gauss.variance,
                        ps=NULL, 
                        Title=Y$Title, 
                        absolute=if (is.null(Y$absolut)) FALSE else Y$absolute,
                        estimate=(if (!is.null(Y$area)) NULL else 
                                  if (is.null(Y$estimate)) TRUE else Y$estimate),
                        ##       currently estimate is always TRUE
                        simulate=(if (!is.null(Y$area)) NULL else
                                  if (is.null(Y$simulate)) TRUE
                                  else Y$simulate),
                        estim.cov=(if (!is.null(Y$area)) NULL else
                                   if (is.null(Y$estim.cov)) simu.standard.estim
                                   else Y$estim.cov),
                        simu.cov=(if (!is.null(Y$area)) NULL else
                                  if (is.null(Y$simu.cov)) simu.standard.cov
                                  else Y$simu.cov), 
                        radius=(if (!is.null(Y$area)) NULL else
                                if (is.null(Y$radius)) simu.radii else Y$radius),
                        scale=if (is.null(Y$scale)) simu.R.vario else Y$scale, 
                        area =Y$area)
          }

          if (!is.null(Y$null.hypo.text)) {
            null.hypo.text <-
              paste("Since $", 
                    paste(simu.name.EVS[Y$null.hypo.text], collapse="$, $"), 
                    "$ is constant, the rejection rates for $", 
                    paste(simu.name.EVS[Y$null.hypo.text], collapse="$, $") , 
                    "$ should be small. ")
          } else  null.hypo.text <- " "

          x.table(matrix(m[[v]][Y$EVS, p, 2:(.mpp.tests - 2)],ncol=.mpp.weights),
                  label=paste(Y$tex.label, Y$label, sep="."), 
                  caption=paste(Y$Title, 
                    if (is.null(Y$no.further.text) || !Y$no.further.text)
                    paste(": average of the rejection rates in \\%", 
                          ##" for the test", 
                          ##ifelse(EVSgiven, paste(" ``$", 
                          ## simu.name.EVS[Y$EVS], "(r)=0$''"), "s"), 
                          " at the ", 
                          100 - Pvalues[p], 
                          "\\% level; the rank of a test among all the ", 
                          length(norm.select) * length(weight.select), 
                          " tests is given in parentheses. ", 
                          sep=""), 
                    null.hypo.text, 
                    Y$caption, sep=""), 
                  file=tex.file.name, 
                  reverse = Y$null.hypo, 
                  norm.select=norm.select, 
                  weight.select=weight.select, 
                  start = (v==1), 
                  end = (v==length(Y$variable)), 
                  add.col=Y$add.col
                  )
          if (is.numeric(dev)) rl("press return")
        }
        if (!file.exists(file.name)) save(file=file.name, m)
      }
    }
  } # latex.teab.123

  latex <- function() {
    ## library(MarkedPointProcess); texfilename <- "srd.jrssb.04"

    stopifnot(.Platform$OS.type=="unix")
    stopifnot(!system("which latex"), !system("which xdvi"))
    texfile <- paste(texfilename, ".tex", sep="")  
    dvifile <- paste(texfilename, ".dvi", sep="")
    f <- scan(paste(systempath, '/tex/', texfile, sep=""),
              what=character(0), sep="\n", blank.lines.skip=FALSE,
              allowEscapes=FALSE)

    print(f)
    
    f[1] <- paste("\\def\\path{", ps.path, "}", sep="")
    f[2] <- paste("\\def\\texpath{", tex.path, "}", sep="")
    f[3] <- paste("\\def\\additive{", simu.models$add$file, "}", sep="")
    f[4] <- paste("\\def\\variance{", simu.models$var$file, "}", sep="")
    f[5] <- paste("\\def\\nn{", simu.models$nn$file, "}", sep="")
    if (f[6]!="") stop("damaged tex file")
    if (file.exists(texfile))
      warning(paste("LaTeX file '", texfile, "' has be overwritten", sep=""))
    write(file=texfile, f)
    file.remove(dvifile)
    system(paste("latex", texfile))
    system(paste("latex", texfile))
    system(paste("xdvi", dvifile))
  }

  #########################################################################
  #########################################################################
  
  items <-
    c("simulate coin model (Ex. 1), intensity=50", 
      "simulate variance model (Ex. 3), intensity=50", 
      "simulate coin model (Ex. 1), intensity=100", #3
      "simulate variance model (Ex. 3), intensity=100", 
      "simulate coin model (Ex. 1), intensity=200", 
      "simulate variance model (Ex. 3), intensity=200", #6
      "simulate nearest neighbour model (Ex. 2)",
      "plot figure 1 (simulations above must have been performed)",
      "tables 1-3 (simulations above must have been performed)", #9
      "test and plot data (Coulissenhieb)", 
      "test and plot data (Steigerwald)", 
      "test and plot data (Coulissenhieb and Biondi et al.)", #12
      "sketches of the theoretical examples",
      "create dvi (#8, #9, #12 must have been performed -- linux machines only!)"
      )

  while (TRUE) {
   if (length(input)==0) input <- menu(items)
    if (input[1] > length(items) || input[1]==0) break
    if (PrintLevel>1) cat(input[1], ":", items[input[1]], "\n")

    switch(input[1], 
           {
             if (!(50 %in% simu.individ.list))
               cat("not considered in the summaries, hence jumped\n")
             else worksheet.add(i=50, w=repet)
           }, {
             if (!(50 %in% simu.individ.list))
               cat("not considered in the summaries, hence jumped\n")
             else worksheet.add(i=50, model=simu.models$var, w=repet)
           }, { #3
             worksheet.add(i=100, w=repet)
           }, {
             worksheet.add(i=100, model=simu.models$var, w=repet)
           }, {
             if (!(50 %in% simu.individ.list))
               cat("not considered in the summaries, hence jumped\n")
             else worksheet.add(i=200, w=repet)
           }, { #6
             if (!(50 %in% simu.individ.list))
               cat("not considered in the summaries, hence jumped\n")
             else worksheet.add(i=200, model=simu.models$var, w=repet)
           }, {
             worksheet.nn(w=repet)
           }, {
             plotFig1(dev=dev)
           }, { #9
             latex.tab.123(dev)
           }, {
             plotdata(data=data[2]) # coulissenhieb only
           }, {
             plotdata(data=data[1]) # steigerwald only
           }, { #12
             if (is.null(biondi.etal)) {
               cat("The data by Biondi et al. are not included in the package.",
                   "\nPlease contact the authors directly and pass the data",
                   "\nto the function in the same format as coulissenhieb",
                   "\nusing the variable biondi.etal. See manual for further",
                   "\ndetails.\n")
             } else plotdata(data=data[c(2,3)]) # biondi only
           }, {
             if (!require(adapt)) next
             theoretical.examples()
           }, {
             try(latex())     
           })
    input <- input[-1]
  } # while (true)
  on.exit(NULL)
  invisible(NULL)
}

