\name{rfm.test}
\alias{rfm.test}
\title{MC test on random field model}
\description{
  \command{rfm.test} performs MC tests which enables the user to decide
  whether a marked point process may be considered as a random field
  model, i.e., as a model where the marks are independent of the locations
%  (but the marks themselves may be dependent, in
%  contrast to the assumption in the Barnard-Besag-Diggle test)
}
\usage{
rfm.test(coord=NULL, data, normalize=TRUE, MCrepetitions=99,  
         MCmodel=list(model="exponential",
           param=c(mean=0,variance=NA,nugget=0,scale=NA)),
         method=NULL,
         bin=c(-1,seq(0,1.2,l=15)), Ebin=seq(0,1,0.01),
         MCregister=1, n.hypo=1000,
         pvalue=c(90, 95, 99), tests="l1 & w3",
         tests.lp=NULL, tests.weight=NULL, Barnard=FALSE,
         PrintLevel=RFparameters()$Print,...
         )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{coord}{matrix with 2 columns; the coordinates of the pointss}
  \item{data}{vector or matrix; the univariate marks that correspond to
      the locations; if \code{data} is a matrix then each column is
      interpreted as an independent observation given the locations
      \code{coord}; see Details for further possibilities
    }
  \item{normalize}{logical; if \code{TRUE} the data are transformed to
    standard normal data before analysed; if \code{data} is a matrix
    this is done for each column separately}
  \item{MCrepetitions}{usually 19 or 99; number of simulations that are
    compared with the data}
  \item{MCmodel}{
    variogram model to be fitted, see
    \command{\link[RandomFields]{fitvario}}.
  }
  \item{method}{method used to simulate Gaussian random fieldsl;
    see \command{\link[RandomFields]{GaussRF}}
  }
  \item{bin}{sequence of increasing bin margins for calculating
    the function E, V, etc; see Details}
  \item{Ebin}{sequence of increasing bin margins for the resulting
    relative MC test positions of the data; see Details}
  \item{MCregister}{0:9; the register to which intermediate results are
    stored when the random fields are generated for the MC test}
  \item{n.hypo}{number of repeated MC tests, see Details}
  \item{pvalue}{test levels}
  \item{tests}{vector of characters, see Details.}
  \item{tests.lp}{vector of characters, see Details.}
  \item{tests.weight}{vector of characters, see Details.}
  \item{Barnard}{test by Barnard (1963) on the independence of marks}
  \item{PrintLevel}{If zero}
  \item{...}{any parameter for \command{\link[RandomFields]{variofit}}
    can be passed,
    except for \code{x}, \code{y}, \code{z}, \code{T}, \code{data},
    \code{model}, \code{param}, \code{mle.methods} and
    \code{cross.methods}
  }
}
\details{
  \code{data}: there are three possibilities to pass the data
  \itemize{
    \item \code{data} a vector or matrix, \code{coord} contains the
    coordinates, as described above

    \item \code{data=list(coord=,data=)} and \code{coord==NULL}

    \item \code{data=list( list(coord=,data=), \ldots, list(coord=,
      data=))}; several data sets are analysed and all the results are
    summed up, and returned in a single matrix \code{E} (or \code{V} or
    \code{SQ})
    }

   \code{bin}: as the variogram in geostatistics, the characteristics for
  the marks of a marked point process depend on a distance (vector)
  \eqn{r}.  Instead of returning a cloud of values, binned values are
  calculated in the same way the binned variogram is obtained.  \code{bin} gives
  the margins of the bins (left open, right closed ones) as an
  increasing sequence.  The first bin must include the zero, i.e.,
  \code{bin=c(-1, 0, \ldots)}.
  \cr
  
  \code{Ebin} is ignored if only a single realisation of the data is given.
  Otherwise \code{Ebin} gives the bounds of the bins for the calculated
  test statistics.

  \code{n.hypo} : the testing algorithm for a data set is as follows:
  \itemize{
    \item the parameters of a Gaussian random field are adapted to the
    data set
    \item \code{n.hypo} realisations of the Gaussian random field are
    are simulated
    \item for each realisation, an 
    MC test based on \code{MCrepetitions} is performed (estimation of
    the parameters of the random field, and test statistics for 99
    realisations)
  }

  \code{tests}, \code{tests.lp}, \code{tests.weight}:
  \itemize{
    \item
    If \code{tests="all"} then the results for all test variants are
    returned, independently of the values of \code{tests.lp} and
    \code{tests.weight}
    
    \item else 
    the results of all tests given by \code{tests} and
    the combinations of \code{tests.lp} and \code{tests.weight}
    are given.

     \item
    The values for \code{tests.lp} are
    \dQuote{max} (maximum norm),
    \dQuote{l2} (l2 norm),
    \dQuote{l1} (l1 norm),
    \dQuote{robust} (the distance is squared for small distances only),
    \dQuote{anti} (the distance is square for large distances only)   
    
    \item
    The values for \code{tests.weight} are
    \dQuote{const} (constant weight),
    \dQuote{1/sum\#} (\sQuote{sum\#} is the cummulative sum of the number of points
    in all bins to the left, and the considered bin itself),
    \dQuote{sqrt(1/sum\#)} (sqrt of \sQuote{1/sum\#}),
    \dQuote{1/sumsqrt\#} (similar to \sQuote{1/sum\#}, but square root of
    the number of points is summed up),
    \dQuote{\#} (number of points within a bin),
    \dQuote{sqrt\#} (square root of the number of points),
    \dQuote{1/sd} (sd=estimated standard deviation within a bin)
    
    or, equivalently,
    
    \dQuote{w1}, \dQuote{w2}, \dQuote{w3}, \dQuote{w4}, \dQuote{w5},
    \dQuote{w6}, \dQuote{w7}.
    
   \item The values for \code{tests} are
    \dQuote{max & const}, \dQuote{l2 & const}, \dQuote{l1 & const},
    \dQuote{robust & const}, \dQuote{anti & const},
    \dQuote{max & 1/sum\#}, \dQuote{l2 & 1/sum\#}, \dQuote{l1 & 1/sum\#},
    \dQuote{robust & 1/sum\#}, \dQuote{anti & 1/sum\#},
    \dQuote{max & sqrt(1/sum\#)}, \dQuote{l2 & sqrt(1/sum\#)},
    \dQuote{l1 & sqrt(1/sum\#)}, \dQuote{robust & sqrt(1/sum\#)},
    \dQuote{anti & sqrt(1/sum\#)}, \dQuote{max & 1/sumsqrt\#},
    \dQuote{l2 & 1/sumsqrt\#}, \dQuote{l1 & 1/sumsqrt\#},
    \dQuote{robust & 1/sumsqrt\#}, \dQuote{anti & 1/sumsqrt\#},
    \dQuote{max & \#}, \dQuote{l2 & \#}, \dQuote{l1 & \#},
    \dQuote{robust & \#}, \dQuote{anti & \#}, \dQuote{max & sqrt\#},
    \dQuote{l2 & sqrt\#}, \dQuote{l1 & sqrt\#}, \dQuote{robust & sqrt\#},
    \dQuote{anti & sqrt\#}, \dQuote{max & 1/sd}, \dQuote{l2 & 1/sd},
    \dQuote{l1 & 1/sd}, \dQuote{robust & 1/sd}, \dQuote{anti & 1/sd},
    
    or, equivalently,

    \dQuote{max & w1}, \dQuote{l2 & w1}, \dQuote{l1 & w1},
    \dQuote{robust & w1}, \dQuote{anti & w1}, \dQuote{max & w2},
    \dQuote{l2 & w2}, \dQuote{l1 & w2}, \dQuote{robust & w2},
    \dQuote{anti & w2}, \dQuote{max & w3}, \dQuote{l2 & w3},
    \dQuote{l1 & w3}, \dQuote{robust & w3}, \dQuote{anti & w3},
    \dQuote{max & w4}, \dQuote{l2 & w4}, \dQuote{l1 & w4},
    \dQuote{robust & w4}, \dQuote{anti & w4}, \dQuote{max & w5},
    \dQuote{l2 & w5}, \dQuote{l1 & w5}, \dQuote{robust & w5},
    \dQuote{anti & w5}, \dQuote{max & w6}, \dQuote{l2 & w6},
    \dQuote{l1 & w6}, \dQuote{robust & w6}, \dQuote{anti & w6},
    \dQuote{max & w7}, \dQuote{l2 & w7}, \dQuote{l1 & w7},
    \dQuote{robust & w7}, \dQuote{anti & w7}

    and
    \dQuote{range} (difference largest positive and largest negative
    deviation for all bins),
    \dQuote{no.bin.sq} (l2 norm where the bins are chosen so that they
    contain only 1 point), 
    \dQuote{no.bin.abs} (l1 norm where the bins are chosen so that they
    contain only 1 point)
  }
}
\value{
  \code{list(E=,VAR=,SQ=,M=,est=,\ldots)} where \code{\ldots}
  are the input parameters such as
  \code{normalize, MCrepetitions, MCmodel, MCparam,
    sill, bin, Ebin}. 
  Let \eqn{n} be the number of currently implemented versions of the MC
  test (using different weights and lp-norms).
  Then \code{VAR}, \code{SQ}, and \code{M} are all matrices with
  \eqn{n} columns. The number of rows depends on the input parameters:
  If only one realisation of the \code{data} is given then
  the absolute test positions of the MC test is returned in \code{E},
  \code{VAR}, \code{SQ}, and \code{M} in a single row.
  If several realisations of the \code{data} (and the \code{coord})
  are given, then the number of rows equals \code{length(Ebin)-1},
  and the each entry contains the number of statistics falling into
  respective (relative) bin given by \code{Ebin}.
}
\references{
  Barnard, G. (1963)
  Discussion paper to M.S. Barlett on \dQuote{The spectral analysis of
    point processes}, \emph{J. R. Statist. Soc. Ser. B}, \bold{25}, 294.

  Besag, J. and Diggle, P. (1977)
  Simple Monte Carlo tests for spatial pattern.
  \emph{J. R. Statist. Soc. Ser. C}, \bold{26}, 327--333.
  
  Schlather, M., Ribeiro, P. and Diggle, P. (2004)
  Detecting Dependence Between Marks and Locations of Marked Point Processes
  \emph{J. R. Statist. Soc., Ser. B} \bold{66}, 79-83.
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\seealso{\command{\link{mpp.characteristics}}, \command{\link{simulateMPP}} }

\examples{}
\keyword{spatial}% at least one, from doc/KEYWORDS

