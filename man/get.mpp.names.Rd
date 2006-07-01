\name{get.mpp.names}
\alias{get.mpp.names}
\title{Model names for marked point process}
\description{
  \code{get.mpp.names} returns the names of implemented marked point processes
}
\usage{
get.mpp.names()
}
\value{
  \command{get.mpp.names} returns a vector of names.
}
\details{
  currently implemented models are
  \itemize{
    \item nearest neighbour\cr
    the points are given by a stationary Poisson point process, and a
    mark is the distance to the nearest neighbour within the process;
    \code{nearest neighbour} has one parameters, multiplied to the
    calculated distance.
    \item random coins\cr
    We start in \eqn{R^d}, here \eqn{d=2},
    with a marked Poisson process where the points are given by
    a stationary Poisson point process \eqn{\Phi}{Phi}
    and the marks are i.i.d.
    random objects (disks or cones) of dimension \eqn{d+1}. 
    At a point \eqn{x} of \eqn{\Phi}{Phi} the sum of the heights of the
    objects that cover \eqn{x} is the mark of this model;
    \code{random coins} has three parameters, the first parameter
    chooses the kind of coin function, the second the scale, the
    third the height.  Currently the available coin functions are
    disk (1) and and cone (2).
    \item variance by coins\cr
    The marks are independent Gaussian random variables with mean 0
    and standard deviation equal to
    the mark of random coin model.
  }
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\seealso{ \command{\link{simulateMPP}} }

\keyword{spatial}
