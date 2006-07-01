\name{srd}
\alias{srd}
\alias{srd.jrssb}
\title{Simulation study for marked point processes}
\description{
  The function reproduces the simulation study published by
  Schlather, Ribeiro and Diggle (2004)
}
\usage{
srd.jrssb(input=NULL, repet=500, dev=2, PrintLevel=2, readlines=TRUE,
          ps.path="ps/", data.save.path="data/", simu.path="simu/",
          tex.path="tex/", biondi.etal=NULL, final=TRUE)
}

\arguments{
  \item{input}{vector of values for automatic choice of the submenues}
  \item{repet}{number of realisations investigated for each parameter
    set. 2 * \code{repet} gives also roughly the number of hours needed
    to do all the simulations!
  }
  \item{dev}{device, see \command{\link[RandomFields]{Dev}}}
  \item{PrintLevel}{If \eqn{\code{PrintLevel}>1} some information on
    the course of the algorithm is printed
  }
  \item{readlines}{logical. Only used if \code{is.numeric(dev)}. If
    the operating system is unix and \code{readlines=FALSE} then the system
    waits a second before the next plot is calculated or shown.
    Otherwise the system waits for return.
  }
  \item{ps.path}{path used for all postscript files in \command{srd.jrssb}
  }
  \item{data.save.path}{directory where intermediate results in the
    analysis of the data are stored}
  \item{simu.path}{directory used to store intermediate results for
    the simulations
  }
  \item{tex.path}{path used for all created tex files containing tables}
  \item{biondi.etal}{
    list of two components:
    \itemize{
      \item{coord:}{a matrix of two columns for
      the coordinates and}
      \item{diameter:}{a vector for the marks}
    }
    The data set by Biondi et al., used in
    Schlather, Ribeiro and Diggle (2004),
    is not freely available.
  }
  \item{final}{logical. If \code{FALSE} then several parameters are changed
    so that the algorithms run faster.
    Especially \code{repet} is limited to 3, \command{\link[stats]{optim}}
    is very sloppy, and the parameter sets for the simulation study
    are remarkably reduced. \code{final=FALSE} is used for debugging
    or demonstration, but the results are not reliable.
  }
}
\details{
  The simulation study is written in a way that \command{srd.jrssb()}
  might be called several times on a parallel CPU system or PVD from the
  same directory.
  The non-interaction of the different instances is managed by means
  of lock files. 
  It is ensured that no simulation study for a specific
  parameter set is performed twice. However,
  the programming has been kept simple so that
  a lock files might be created but the simulation is
  not performed by any process because of strange
  interaction. In this case all remaining lock files have
  to be deleted after the first run of the simulation study, and
  the study has to be re-run. Then, \command{srd.jrssb()} does the simulation
  study only for remaining cases.

  For both the simulation studies and the data analyses, intermediate
  stored results are re-used if run again.
}

\value{
  \code{NULL}.
  Side effect is the creation of postscript files and tex files.
}
\references{
  Schlather, M., Ribeiro, P. and Diggle, P. (2004)
  Detecting Dependence Between Marks and Locations of Marked Point Processes.
  \emph{J. R. Statist. Soc., Ser. B} \bold{66}, 79-83.

  Biondi, F., Myers, D.E., and Avery, C.C. (1994)
  Geostatistically modeling stem size and increment in an old-growth
  forest.
  \emph{Can. J. Forest Res.} \bold{24}, 1354-1368.
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}

\seealso{ \command{\link{rfm.test}}, \command{\link{simulateMPP}},
  \code{\link{MarkedPointProcess}}
}

\examples{
% library(MarkedPointProcess, lib="~/TMP"); 
# repet=2 is for demonstration of functionality
# repet=50 gives more reliable results
# repet=500 and final=TRUE have been used in the simulation study
  
printlevel <- 1 + !interactive()
srd.jrssb(input=if(!interactive()) c(1:2,0), repet=2, dev=TRUE,
          final=FALSE, Pr=printlevel)
}
\keyword{spatial}


% library(RandomFields, lib="~/TMP"); Dev(TRUE,TRUE,ps="xx");plot(1,1);Dev(FALSE)
