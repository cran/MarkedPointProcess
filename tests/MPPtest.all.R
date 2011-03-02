## R --no-save < MPPtest.all.R
# source("MPPtest.all.R")

if (EXTENDED.TESTING <- file.exists(f <- "MPPtest.R")) source(f) else 
if (file.exists(f <-"~/R/MPP/MarkedPointProcess/tests/MPPtest.R")) source(f) else
q()

RFparameters(PrintLevel=1, PracticalRange=FALSE)

repetitions <- 1
radius <- 0.5
range <-0.5

if (EXTENDED.TESTING) {
  coordrepetitions <- 2000
  gaussvariance <-c(0.000001,(1:10)/10);# gaussvariance <- 1.0
  radius <- c(0.01,0.05,0.1,0.5)
  range <-  c(0.01,0.05,0.1,0.5,1,2)
  coordmodels <- c("uniform", "Poisson")
} else {
  coordrepetitions <- 20
  gaussvariance <-c(0.000001,(1:2)/2);# gaussvariance <- 1.0
  radius <- c(0.01,0.1)
  range <-  c(0.05,0.1)
  coordmodels <- c("Poisson")
}

model <- 'list("+", list("$", var=gv, scale=gs, list("expon")),\
               list("random coins", param=c(1, ar, ah)))'
compare.model <- 'list("+", list("$", var=gv, scale=gs, list("expon")), \
                       list("$", var=(sill - gv), scale=2 * ar, list("circ")))'


unit <- 2
individsperunit2 <- c(50)

xlim <- c(0,1) * unit
ylim <- c(0,1) * unit
edgecorrection <- TRUE; # edgecorrection <- FALSE
endofbins <- 0.7


for (individ in individsperunit2) {
  for (coordmodel in coordmodels) {
    x <- MPPcontrol(individ=individ * unit * unit,
                    lambda = individ,
                    repetitions=repetitions,
                    edgecorrection = edgecorrection,
                    coordmodel=coordmodel,coordrepetitions=coordrepetitions, 
                    gauss.variance=gaussvariance, gauss.scale=range,
                    rc.radius=radius, #nn.factor=1,
                    model = model,
                    compare.model=compare.model,
                    endofbins=endofbins,numberbins=30, 
                    xlim=xlim,ylim=ylim, hist=FALSE, wait=FALSE,
                    waitfinally=FALSE)
  }
}


