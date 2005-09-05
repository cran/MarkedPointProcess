
if (file.exists("source.R")) source("source.R")

##compare.parameter="c(0,gv+gs^2,ar+ah,sqrt(nf))",
MPPcontrol <- function (individ=50, lambda =50, coord=NULL,
                        coordmodel=c("uniform", "given", "Poisson"),
                        repetitions=100, coordrepetitions=1,
                        edgecorrection=FALSE,
                        gauss.variance=NA, gauss.scale=NA,
                        rc.radius=NA, nn.factor=NA,    
                        model,
                        sill = 1,
                        compare.model=NULL,
                        endofbins=0.7, numberbins=30,
                        xlim=c(0,1), ylim=c(0,1),
                        hist=FALSE, wait=FALSE, clean=FALSE, waitfinally=TRUE
                       ){
  coordmodel = match.arg(coordmodel)
  ## for general structure see RFcontrol
  ##
  ## model : string which can be evaluated to list(list(),op,...list())
  ##         containing the parameters gv(gauss variance), gs(scaling),
  ##                       ar(radius of random coins), ah(calculated height),
  ##                       nf(nn.factor) 
  ##
  ## gauss.variance, gauss.scale, rc.radius, nn.factor might be
  ##         vectors; then it is looped through all combinations
  ##
  ## compare.model: if NULL no comparision model else compare.model is of
  ##                string that can be evaluated to list(list(),op,...,list())
  dim <- 2
  digits <- 2
  coordmodel <- match.arg(coordmodel)
  
  if (is.null(compare.model)) {
    compare.variogram <- function(x, gv, gs, ar, nf) { return(NA); }
  } else {
    compare.variogram <- function(x, gv, gs, ar, nf) {
      ## that time `Variogram' was not programmed yet
      CovarianceFct(0, model=eval(parse(text=compare.model)), dim=dim) -
      CovarianceFct(x, model=eval(parse(text=compare.model)), dim=dim)
    }
  }

  if (coordmodel=="Poisson")
    { npoints <- NULL } else
  if (coordmodel=="uniform") {npoints <- individ} else
  stop("<given> not allowed here")
   
  modelname <- paste(as.character(model[[1]]),collapse=",")
  if (length(model)>1)
    for (i in 2:length(model))
      modelname <-
        paste(modelname, "_", paste(as.character(model[[1]]),collapse=","))

  meanvar <- "" ## to avoid an error in the final output in case everything fails
  for (gv in gauss.variance)  for (gs in gauss.scale)
    for (ar in rc.radius) for (nf in nn.factor){    
      edgecor <- ar
      ah <- -sqrt((sill-gv) * diff(xlim)*diff(ylim)/(individ*pi*ar^2))
      rc.parameter <- c(rc.radius=ar,rc.height=ah)
      nn.parameter<-c(nf)
      
      modeletc <-paste(modelname,
                       ",gv=",format(gv,dig=digits),
                       ",gs=",format(gs,dig=digits),
                       ",ar=",format(ar,dig=digits),
                       ",ah=",format(ah,dig=digits),
                       ",nf=",format(nf,dig=digits),
                       ",ind=",individ,
                       ",co.me=",coordmodel,
                       sep="")
      print(modeletc)
  
      bin <- c(-1,seq(0,endofbins,l=numberbins))  
      midbin <- as.double(0.5 * (bin[-length(bin)] + bin[-1]));
      midbin[1] <- 0
      truecov <- double(numberbins)
      binresult <- double(numberbins)
    
      v <- seq(0,0,l=numberbins)
      repet <- seq(0,0,l=numberbins)
      print("start simu")
      simu <- simulateMPP(coord=coord, npoints=npoints, repetitions=repetitions,
                           coordmodel=coordmodel, window=c(xlim, ylim),
                           edgecorrection=if (edgecorrection) edgecor else 0,
                           lambda=lambda,
                           coordrep=coordrepetitions,
                           model=eval(parse(text=model)))
      print("end simu")
      ##print(simu)
      if (is.numeric(simu)) {print(c("simu return value:",simu)); break;}
    
      for (i in 1:coordrepetitions) {
        ##print(c("coordrepet",i))
        allres <- simu[[i]][[2]]
        meanvar <-
          paste("m=", format(mean(allres), dig=digits),
                ", sd=", format(sqrt(var(as.vector(allres))),dig=digits))
        if (hist) hist(as.vector(allres))
        binresult <-
          EmpiricalVariogram(x=simu[[i]]$coord[,1], y=simu[[i]]$coord[,2],
                             data=allres, grid=FALSE, bin=bin)$emp.var
        index <- is.finite(binresult) & binresult >= 0;
        v[index] <- v[index] + binresult[index];
        repet[index] <- repet[index] + 1;
      } # coordrep
      v <- v / repet;
      truevario <- compare.variogram(midbin, gv, gs, ar, nf)
      ##print(truevario);
      delta <- mean(abs(v-truevario),na.rm=TRUE)
      print(delta)
    
      if (!is.na(delta)) plot(midbin,v,type="l",ylim=c(0,max(1,max(truevario))))
      else plot(midbin, -1, type="l", ylim=c(0, max(1, max(truevario))))
      points(midbin,truevario)
      if (!interactive())
        title(sub=paste("delta=", paste(format(delta, dig=digits), collapse=","),
                ",",meanvar,sep=""),
              main=paste(modeletc))
      if (wait) readline()
    }
  if (waitfinally && !wait)  { print("press key"); readline(); }
  if (clean) while (!is.null(dev.list())) dev.off();
  return(NULL);
}

