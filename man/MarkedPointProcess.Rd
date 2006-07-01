\name{MarkedPointProcess}
\alias{MarkedPointProcess}
\title{Simulation study for marked point processes}
\description{
  This package allows for simulating and analysing marked point processes
 }

 \details{
   The following functionalities are provided:
   \itemize{
     \item \command{\link{get.mpp.names}} : returns the names of the
     implemented models
     \item \command{\link{mpp.characteristics}} : returns 
     characteristics for the marks of marked point processes such as
     the mark variogram, Stoyan's kmm function, the E function and the
     V function
     \item \command{\link{rfm.test}} : MC test whether E or V is a constant.
     If any of these hypotheses are rejected, the investigated marked
     point process cannot be considered as random field model, i.e. a
     model where the marks are independent of the locations (however the
     random field model allows that the marks themselves are spatially
     dependent) 
     \item \command{\link{simulateMPP}} : simulation of marked point
     processes
     \item \command{\link{srd.jrssb}} : function that generates the results
     published by Schlather, Ribeiro, Diggle (2004)
     \item \command{\link{splitmodel}} : auxilliary function that splits a
     user defined model in a pure Gaussian random field part and a pure
     marked point process part
   }
   Further, a forestry data set is provided, see \link{BITOEK}.
}

\references{
%  \itemize{
%    \item
    Schlather, M., Ribeiro, P. and Diggle, P. (2004)
    Detecting Dependence Between Marks and Locations of Marked Point Processes
    \emph{J. R. Statist. Soc., Ser. B} \bold{66}, 79-83.
%    \item
%    Stoyan, D., Kendall, W.S. and Mecke, J. (1995)
%    \emph{Stochastic Geometry and its Applications.}
%    Chichester: Wiley.
}

\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}

\section{Acknowledgement}{
  The work has been financially
  supported by the German Federal Ministry of Research
  and Technology (BMFT) grant PT BEO 51-0339476C during 2000-03.
}

\keyword{spatial}


