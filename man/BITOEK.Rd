\name{BITOEKforests}
\alias{BITOEK}
\alias{BITOEKforests}
\alias{forests}
\alias{steigerwald}
\alias{coulissenhieb}
\title{Forestry areas of intensif measurement by BITOEK}
\usage{data(BITOEK)}
\description{
  Two foresty areas of intensive measurement maintained by BITOEK
}
\format{ 
  The two variables, \code{steigerwald} and \code{coulissenhieb}
  are lists, where the first component, \code{coord},
  is a matrix of two columns giving
  the positions of the trees, and the second component gives the
  stem diameters. The third component, \code{author},
  gives a short description of the data set.
}
\describe{
  Both forests can be found in the the North of Bavaria, Germany;
  \sQuote{Steigerwald} is a forest near Bamberg, \sQuote{Coulissenhieb}
  is in the \sQuote{Fichtelgebirge}.
}
%\details{
%}
\source{
  The data were collected by Pedro Gerstberger,
  BITOEK, \url{http://www.bitoek.uni-bayreuth.de},
  University of Bayreuth, Germany.
}
\references{
  Schlather, M., Ribeiro, P. and Diggle, P. (2004)
  Detecting Dependence Between Marks and Locations of Marked Point Processes
  \emph{J. R. Statist. Soc., Ser. B} \bold{66}, 79-83.
}
\examples{
%  library(MarkedPointProcess, lib="~/TMP")
data(BITOEK)

plotWithCircles(cbind(coulissenhieb$coord, coulissenhieb$diam), factor=3)

plotWithCircles(cbind(steigerwald$coord, steigerwald$diam), factor=2)

}
\keyword{datasets}






