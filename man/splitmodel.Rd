\name{splitmodel}
\alias{splitmodel}
\title{Split betweem marked point processes and random fields }
\description{
  \command{splitmodel} splits a model given in form of a list
  (the third variant of model definition for random fields, see
  \command{\link[RandomFields]{CovarianceFct}}) into a random field
  part and a marked
  point process part
}
\usage{
splitmodel(model)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{model}{The definition of a model is of the form
    \code{model = list(l.1, OP.1, l.2, OP.2, \ldots, l.n)}.
    The lists \code{l.i}
    are all either of the  form \code{l.i = list(model=,var=,kappas=,scale=)}
    or of the form \code{l.i = list(model=,var=,kappas=,aniso=)} in case
    of random field parts, or of the form \code{l.i =
      list(model=,param=)} in case of marked point process parts.
    \code{l.i$model} is a string; \code{var} gives the variance;
    \code{scale} is a scalar whereas \code{aniso} is a \eqn{d \times
      d}{d x d} matrix, which is multiplied from the right to the points, and
    at the transformed points the values of the (isotropic) random field
    (with scale 1) are
    calculated.  The dimension \eqn{d} of matrix must match the
    number of rows of \code{x}.  \code{param} is vector of real values
    whose length depends on the specified \code{model}.  The
    models for the random field part
    can be combined by \code{OP.i="+"} or \code{OP.i="*"}, those for the
    marked point process parts only by \code{OP.i="+"}.
  }
}
\value{
  \code{list(RF=RF, mpp=mpp)} where \code{RF} is a usual model
  definition for a random field.  Further,
  \code{mpp=list(mpp.1,\ldots,mpp.n)},
  where \code{mpp.i=list(model=model,param=param,mnr=)} and \code{mnr}
  is the internal C code for \code{model}.
}

\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}
\seealso{\command{\link{simulateMPP}} }

\examples{
str(splitmodel(list(list(model="exp", var=5, scale=3))))

str(splitmodel(list(list(model="nearest neighbour", param=4))))

str(splitmodel(list(list(model="exp", var=5, scale=3),
                    "+",
                    list(model="nearest neighbour", param=4)
                    )))

str(splitmodel(list(list(model="exp", var=5, scale=3),
                    "*",
                    list(model="spherical", var=1, scale=2),
                    "+",
                    list(model="nearest neighbour", param=4),
                    "+",
                    list(model="random coin",
                         param=c(fct=1, scale=7, height=8))
                    )))
}
\keyword{spatial}% at least one, from doc/KEYWORDS

