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
         bin=c(-1,seq(0,1.2,l=15)), MCregister=1, n.hypo=1000,
         pvalue=c(10, 5, 1), tests="l1 & w3",
         tests.lp=NULL, tests.weight=NULL, Barnard=FALSE,
         PrintLevel=RFparameters()$Print,...
         )
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{coord}{matrix with 2 columns; the coordinates of the points}
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
  \item{method}{method used to simulate Gaussian random fields;
    see \command{\link[RandomFields]{GaussRF}}
  }
  \item{bin}{sequence of increasing bin margins for calculating
    the function E, V, etc in analogy to the binning for variograms;
    see Details}
  \item{MCregister}{0:9; the register to which intermediate results are
    stored when the random fields are generated for the MC test}
  \item{n.hypo}{number of repeated MC tests to determine the
    \code{pvalue}-position for the Nullhypothesis.
    If the variogram were not estimated, this position would
    be \eqn{(1-\code{pvalue}) (\code{MCrepetitions} + 1)}.
    see Details}
  \item{pvalue}{test levels in percent.
    Only values below 50 are accepted;
    otherwise 100-\code{pvalue} is regarded as \code{pvalue}
    (to be consistent with the former definition)
  }
  \item{tests}{vector of characters, see Details.}
  \item{tests.lp}{vector of characters, see Details.}
  \item{tests.weight}{vector of characters, see Details.}
  \item{Barnard}{test by Barnard (1963) on the independence of marks}
  \item{PrintLevel}{If zero then no messages are printed. The higher the
  value the more tracing information is given.}
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
 
  \code{n.hypo} : for correct appreciation of the relative position of the
  statistic for the data set w.r.t. the simulations,
  the reference values for the estimated \code{pvalue} level must be determined:
  \itemize{
    \item the parameters of a Gaussian random field are adapted to the
    data set
    \item \code{n.hypo} realisations of the Gaussian random field are
    are simulated
    \item for each realisation, an 
    MC test based on \code{MCrepetitions} is performed (estimation of
    the parameters of the random field, and test statistics for 99
    realisations)
    \item the position of the
    statistic for the simulation within the statistics of
    the 99 realisations are calculated.
    \item
    The position of the larges \code{pvalue} positions is determined, which
    have values around \eqn{(1-\code{pvalue}) (\code{MCrepetitions} +
      1)};
    this position is returned in \code{null.sl} as
    reference values for the estimated \code{pvalue} level
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
  Let \eqn{n} be the number of MC tests chosen by the user.
  Then \command{rfm.test} returns a list of the following elements:
  \cr

  \code{E} (\eqn{E} function)\cr
  matrix of \eqn{n} columns. The number of rows depends on the input
  parameters:
  If only one realisation of the \code{data} is given then
  the absolute test positions of the MC test is returned, i.e. a
  value between 1 and \eqn{\code{MCrepetition}+1}, inclusively.
  
  If several realisations of the \code{data} (and the \code{coord})
  are given, then the number of rows equals \eqn{\code{MCrepetition}+1},
  and the \eqn{k}th row gives
  the number of test statistics with position \eqn{k}.
  The first situation is the standard one for the user. The second
  situation appears when \command{rfm.test} is recalled to calculate
  the intermediate result \code{null.hypo}, see below.
  \cr
     
  \code{VAR} (\eqn{V} function)\cr
  matrix of \eqn{n} columns. See \code{E} above. 
  \cr

  \code{SD} (the square root of the \eqn{V} function)\cr
  matrix of \eqn{n} columns. See \code{E} above. 
  \cr

  \code{reject.null}\cr
  list of logical matrices that indicate whether \code{E}, \code{VAR}
  or \code{SQ} should be rejected at the given levels, i.e. whether
  the positions of the tests statistics for \code{E}, \code{VAR}
  or \code{SQ} are at least as large as the estimated reference values
  given by \code{null.sl}.
  \cr
  
  \code{est}\cr
  list of variogram models according to \code{MCmodel} estimated 
  from the data.
  \cr

  \code{normalize}\cr
  The input parameter \code{normalize}.
  \cr

  \code{MCrepetitions}\cr
  The input parameter \code{MCrepetitions}.
  \cr

  \code{MCmodel}\cr
  The input parameter \code{MCmodel}.
  \cr

  \code{null.hypo}\cr
  \code{null.hypo} stores intermediate results that are usually not of
  interest for the user.
  \code{n.hypo} simulations have been performed under the null
  hypothesis to determine the \code{pvalue} test positions. (The explicite
  determination is necessary, since parameters of the variogram
  have to be estimated within the null hypothesis.)
  For these \code{n.hypo} simulations, \command{rfm.test} is run and
  \code{null.hypo} gives the results. Note that here, all test variants
  are considered.
  \cr
  
  \code{null.sl}\cr
  List of matrices. They  give the
  reference values for the estimated \code{pvalue} level.
  The values are around (1-\code{pvalue})*(\code{MCrepetitions} + 1),
  but can range between 1 and \eqn{\code{MCrepetitions} + 2}.
  If a value of \eqn{\code{MCrepetitions} + 2} occurs, usually,
  \code{MCrepetitions} and/or \code{n.hypo} have been chosen too small.
  \cr
  
  \code{bin}\cr
  the binning used to calculate \code{E}, \code{VAR} and \code{SQ}
}
\note{
   In comparison to version 0.1 of MarkedPointProcess and the paper by
   Schlather et al. (2004), the announced positions of the test statistics
   for \code{E}, \code{VAR},
   \code{SD} and \code{null.sl} are all increased by \eqn{1},
   now ranging from 1 to 100 instead of 0 to 99, for the standard settings.
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

\examples{
%  library(MarkedPointProcess, lib="~/TMP")
data(BITOEK)
d <- steigerwald
plotWithCircles(cbind(d$coord, d$diam), factor=2)
mpp.characteristics(x=d$coord, data=d$diam,
                    bin=c(-1, seq(0, 50, 2)), show=interactive())

## testing for E=const, V=const or SD=const (this takes several minutes!)
res <- rfm.test(d$coord, d$diam, MCrep=if (interactive()) 99 else 9,
                n.hypo=if (interactive()) 100 else 2)
% res <- rfm.test(d$coord, d$diam, MCrep=9, n.hypo=2, tests="all")

## test statistics for the data
res$E
res$VAR

## reference values for the estimated 10\%, 5\% and 1\% level
res$null.sl

## should E=const, V=const or SD=const be rejected at the given levels?
res$reject.null 
}

\keyword{spatial}% at least one, from doc/KEYWORDS

%% library(RandomFields); source("~/R/RF/RandomFields/tests/source.R")
%x <- 1:100
%z <- GaussRF(x, model="exp", param=c(0,2,0,5), n=10)
%v <- fitvario(x=x, data=z, cross.methods=NULL, mle.methods="ml",
%              model=list(model="exp", param=c(0, NA, 0, NA)), table.f=TRUE)
