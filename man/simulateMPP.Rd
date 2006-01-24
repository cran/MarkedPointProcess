\name{simulateMPP}
\alias{simulateMPP}
\title{Simulation of marked point processes}
\description{
  \code{simulateMPP} generates realisations of marked point processes
}
\usage{
simulateMPP(coordmodel=c("given", "uniform", "Poisson"),
             coord=NULL, npoints=NULL, lambda=NULL,
             window=NULL,  edgecorrection=0.0,
             repetitions=1, coordrepet=1, model=NULL,
             register=0, method=NULL) 
}
\arguments{
  \item{coordmodel}{if \code{coordmodel="given"} then \code{coord}
    are expected to be given and not simulated;
    if \code{coordmodel="uniform"} then
    \code{npoints} uniformly distributed points are created; if
    \code{coordmodel="Poisson"} then a conditional Poisson point process is
    simulated with intensity \code{lambda}}

  \item{coord}{matrix with 2 columns; coordinates of the points;
    \code{coord} is given only if \code{coordmodel="given"}}

  \item{npoints}{number of coordinates;
    \code{npoints} must be given if \code{coordmodel="uniform"}.
  }

  \item{lambda}{intensity of the Poisson process;
    \code{lambda} must be given if \code{coordmodel="Poisson"}.}

  
  \item{window}{\code{= c(xlim, ylim)}. \code{window} must be given
    if \code{coordmodel} equals \code{"uniform"} or \code{"Poisson"}.
  }

  \item{edgecorrection}{double. If \eqn{\code{edgecorrection}>0}
    then a Poisson process is simulated
    with intensity \code{lambda} in a frame of thickness \code{edgecorrection}
    around the \code{window}.
    If \code{window} is not given, the range of the x values and the
    range of the y values are taken to define the \code{window}.
    If \code{lambda} is not given, the intensity within the
    \code{window} is used instead.}
  
  \item{repetitions}{integer; number of independent drawings of the
    marks for a given set of coordinates}
  
  \item{coordrepet}{number of independent drawing of the coordinates;
    this parameter is ignored in case of \code{coordmodel="given"}}
  
  \item{model}{list of lists; model for the marks; see Details
    and \command{\link{get.mpp.names}}.
  }
  
  \item{register}{the register where intermediate results in the
    Gaussian random field simulation are stored, see
    \command{\link[RandomFields]{GaussRF}}}
  
  \item{method}{the method by which the Gaussian random field is
    simulated; if \code{is.null(method)} then the method is chosen
    automatically, see \command{\link[RandomFields]{GaussRF}}
  }
}
\details{
  The definition of a model is of the form
  \code{model = list(l.1, OP.1, l.2, OP.2, \ldots, l.n)}.
  The lists \code{l.i}
  are all either of the  form \code{l.i = list(model=,var=,kappas=,scale=)}
  or of the form \code{l.i = list(model=,var=,kappas=,aniso=)} in case
  of random field parts, or of the form \code{l.i =
    list(model=,param=)} in case of marked point process parts.
  \code{l.i$model} is a string; \code{var} gives the variance;
  \code{scale} is a scalar whereas \code{aniso} is a \eqn{d \times
    d}{d x d} matrix, which is multiplied from left to the points, and
  at the transformed points the values of the random field are
  calculated.  The dimension \eqn{d} of matrix must match the
  number of rows of \code{x}.  \code{param} is vector of real values
  whose length depends on the specified \code{model}.  The
  models for the random field part
  can be combined by \code{OP.i="+"} or \code{OP.i="*"}, those for the
  marked point process parts only by \code{OP.i="+"}.
}
\value{
  \code{coordrepet==1}\cr
  the function returns \code{list(coord, data)},
  \code{data} contains the independent drawing of the marks (as
  columns)\cr
  
  \code{coordrepet>1}\cr
  the function returns
  \code{list( list(coord, data), \ldots, list(coord, data) )}
}
%\references{}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\seealso{\command{\link{get.mpp.names}},
  \command{\link{rfm.test}}, \command{\link{simulateMPP}},
  \command{\link{splitmodel}}, \link{MarkedPointProcess}
}

\examples{
 xlim <- c(0, if (interactive()) 200 else 20)
 mpp <- simulateMPP(coordmodel="Poisson", lambda=1,
              window=c(xlim=xlim, ylim=c(20, 70)),
              repet=3, coordrepet=4,
              model=list(list(model="exp", var=1, scale=10),
                         "+",
                         list(model="nearest neighbour", p=1)))
  str(mpp)
}
\keyword{spatial}
