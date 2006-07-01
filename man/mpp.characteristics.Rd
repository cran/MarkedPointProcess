\name{mpp.characteristics}
\alias{mpp.characteristics}
\title{Characteristics of the marks of a marked point process}
\description{
  \command{mpp.characteristics} returns summary statistics for the marks of
  a marked point process
}
\usage{
mpp.characteristics(...,bin=NULL, rep=1, p=0.8, name="", normalize=TRUE,
                    show=FALSE, model=NULL, param=NULL, summarize=TRUE,
                    PrintLevel=RFparameters()$Print,
                    dev=if (name=="") 2 else FALSE,
                    rdline=if (is.logical(dev)) NULL else readline,
                    staticchoice=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{\dots}{coordinates and the data for the species, see Details}
  \item{bin}{sequence of increasing bin margins for the functions
    \eqn{E}, \eqn{V}, and \eqn{S}; see Details}
  \item{rep}{number of independent measurements of the marks (at each point); usually 1}
  \item{p}{in \eqn{[0,1]}; outlier threshold for the robustified
    distance function of the test statistics} 
  \item{name}{character; if \code{show=FALSE} this parameter is
    ignored. Otherwise, 
    if \code{name=""} then plots are printed on the screen else \code{name}
    is the main name of the postscript files}
  \item{normalize}{logical; if \code{TRUE} the data are (marginally)
    transformed to Gaussian variables (for each species, each kind of
    mark, and each realisation, seperately)
    before being analysed}
  \item{show}{logical; if \code{TRUE} the results are also shown
    graphically }
  \item{model}{the variogram model to compare with the mark variogram; see
    \command{\link[RandomFields]{CovarianceFct}}}
  \item{param}{the parameters for the variogram model; see
    \command{\link[RandomFields]{CovarianceFct}}}
  \item{summarize}{logical; if \code{FALSE} results are for each
    realisation seperately (instead of being averaged over the realisations).} 
  \item{PrintLevel}{0,1 or 2. The function gives some short messages if
    \code{PrintLevel} is 1 or 2.}
  \item{rdline}{\code{NULL} or function.
    if not \code{NULL} then the function is called after each plot with
    a string parameter that gives file or image information
  }
  \item{dev}{the graphical device for the output, see
    \command{\link[RandomFields]{Dev}}}
  \item{staticchoice}{
    logical. The calculation of some variances requires
    the splitting of the data into group. If \code{staticchoice=FALSE}
    this is done in a random way. \code{staticchoice=TRUE} is only
    used for internal testing.
  }
}
\details{
  \code{bin}: analogously to
  the variogram in geostatistics, the characteristics \eqn{E}, \eqn{V},
  and \eqn{S} of a stationary and isotropic
  marked point process depend on the distance
  \eqn{r}.  Instead of returning a cloud of values, binned values are
  calculated in the same way the binned variogram is calculated.
  \code{bin} gives
  the margins of the bins (left open, right closed ones) as an
  increasing sequence.  The first bin must include the zero, i.e.,
  \code{bin=c(-1, 0, \ldots)}.\cr

  \code{\ldots}: data for mark point processes typically split up into
  different species (ill/healthy cells; beaches/oaks/pines);
  furthermore,
  multivariate data are measured for each individual (size of the cell;
  diameter of the stem, height of the tree).  The function 
  calculates many cross-statistics; for example the cross variogram
  of mark \eqn{A} of species \eqn{B} and mark \eqn{C} of
  species \eqn{D} given species \eqn{B} and species \eqn{D}
  are a distance \eqn{r} apart.
  
  Denote by \eqn{S_i} species \eqn{i}, \eqn{i=1,\ldots,s}{i=1,...,s}.
  Due to the potential complexity of the data,
  the data are passed to \command{mpp.characteristics} in the following way:\cr
  (coordinates of species \eqn{S_i}), (marks of species A),\ldots,
  (coordinates of species Z), (marks of species Z)\cr
  The coordinates are \eqn{(n_i \times 2)}{n_i x 2} matrices; the data
  are \eqn{(n_i \times m \times \code{rep})}{n_i x m x \code{rep}} matrices,
  if the data are \eqn{m}-variate, and \code{rep} independent observations
  of the data exist.  In case \eqn{\code{rep}>1} the sequence for the
  data is:\cr
  columns 1:m : first set of the \eqn{m}-variate data,\ldots,\cr
  columns (m * \code{rep} - m + 1) : (m * \code{rep}) : last set of the
  \eqn{m}-variate data.

  Note that \eqn{m} and \code{rep} must be identical for all species.

  The function returns the following values if \code{summarize=TRUE}.
  Denote by \eqn{M_k(S_i)} the
  \eqn{k}th mark of species \eqn{S_i}, and by \eqn{d_{ij}}{d_ij}
  the distance of two indviduals of species \eqn{S_i} and \eqn{S_j}.
  Denote by \eqn{s} the number of species.
  \itemize{
    \item
    \code{E} : (\code{length(bin)-1})\eqn{\times}{ x }(\eqn{ms^2}) matrix.
    function \eqn{E} in the following ordering of the
    columns.
    The index for the marks runs fastest, then the index of the
    conditioning species, then the index for the species the marks
    belong to. 
    That is,
     \cr
     \eqn{E(M_1(S_1)|d_{11}=r)}{E[M_1(S_1) | d_11=r]},
     \eqn{E(M_2(S_1)|d_{11}=r)}{E[M_2(S_1) | d_11=r]}, ...,
     \eqn{E(M_m(S_1)|d_{11}=r)}{E[M_m(S_1) | d_11=r]},  ...,
     \cr
     \eqn{E(M_1(S_1)|d_{12}=r)}{E[M_1(S_1) | d_12=r]}, ...,
     \eqn{E(M_m(S_1)|d_{12}=r)}{E[M_m(S_1) | d_12=r]},  ...,
     \cr
     \eqn{E(M_1(S_1)|d_{1s}=r)}{E[M_1(S_1) | d_1s=r]}, ...,
     \eqn{E(M_m(S_1)|d_{1s}=r)}{E[M_m(S_1) | d_1s=r]},
     \cr...,
     \cr
     \eqn{E(M_1(S_2)|d_{21}=r)}{E[M_1(S_2) | d_21=r]}, ...,
     \eqn{E(M_m(S_2)|d_{21}=r)}{E[M_m(S_2) | d_21=r]},  ...,
     \cr
     \eqn{E(M_1(S_2)|d_{2s}=r)}{E[M_1(S_2) | d_2s=r]}, ...,
     \eqn{E(M_m(S_2)|d_{2s}=r)}{E[M_m(S_2) | d_2s=r]},
     \cr..., 
     \cr
     \eqn{E(M_m(S_s)|d_{ss}=r)}{E[M_m(S_s) | d_ss=r]}.
     \item
     \code{ETest} : matrix of 66 rows and the same structure of columns as
     \code{E}.  
     The rows correspond to different algorithms for calculating the
     deviance of \code{E} from a horizontal line.  Let denote
     by \eqn{E(i)} the value of the \eqn{i}th bin, \eqn{i=0,...,b},
     where \eqn{E(0)}
     is the bin that includes 0.
     \cr
     1 : \eqn{\max E(i) - \min E(i)}{max E(i) - min E(i)}\cr
     2-36 : see below\cr
     37 : \eqn{\sum |E_{i,unbinned} - E(0)|^2}{sum |E_{i,unbinned} - E(0)|^2}\cr
     38 : \eqn{\sum |E_{i,unbinned} - E(0)|}{sum |E_{i,unbinned} - E(0)|}\cr
     \cr
     \bold{2-37:}\cr
     results of (5 \dQuote{norms}) x (7 set of weights)\cr
     \bold{\dQuote{norms}}\cr
     1 : \eqn{\max w_i |E(i) - E(0)|}{max w_i |E(i) - E(0)|}\cr
     2 : \eqn{\sum w_i |E(i) - E(0)|^2}{sum w_i |E(i) - E(0)|^2}\cr
     3 : \eqn{\sum w_i |E(i) - E(0)|}{sum w_i |E(i) - E(0)|}\cr
     4 : \eqn{\sum w_i * r(E(i) - E(0))}{sum w_i * r(E(i) - E(0))}\cr
     5 : \eqn{\sum w_i * s(E(i) - E(0))}{sum w_i * s(E(i) - E(0))}\cr
     where \eqn{r} is a function that
     first increases quadratically, then linearly:
     \eqn{r(x) = x^2} if \eqn{x < q}, and \eqn{2 q x - q^2}
     otherwise, where \eqn{q} is the \code{p}th quantile of
     { |E(i) - E(0)| : i=0,...,b }. 
     The functions \eqn{s} first increase linearly,
     then quadratically: \eqn{s(x) = x} if \eqn{x<a}
     and, \eqn{b (x+a)^2} otherwise.
     Here, \eqn{a = 1/10} and \eqn{b_k = 0.25 / a}.
     
     \bold{Weights}\cr
     a : \eqn{w_i} ~ 1\cr
     b : \eqn{w_i} ~ \eqn{1 / \sum_{j=0}^{i-1} \code{Ebin}(j)}{1 /
       sum_{j=0}^{i-1} \code{Ebin}(j)} \cr
     c : \eqn{w_i} ~ \eqn{\sqrt{1 / \sum_{j=0}^{i-1}
	 \code{Ebin}(j)}}{sqrt(1 / sum_{j=0}^{i-1} \code{Ebin}(j))} \cr
     d : \eqn{w_i} ~ \eqn{1 / \sum_{j=0}^{i-1} \sqrt{\code{Ebin}(j)}}{1 /
       sum_{j=0}^{i-1} \sqrt(\code{Ebin}(j))} \cr
     e : \eqn{w_i} ~ \eqn{\code{Ebin}(i)}{\code{Ebin}(i)}\cr
     f : \eqn{w_i} ~ \eqn{\sqrt{\code{Ebin}(i)}}{sqrt(\code{Ebin}(i))}\cr
     g : \eqn{w_i} ~ \eqn{1 / \sqrt{var(E_j,unbinned; distance in ith
	 bin)}}{1 / sqrt(var(E_j,unbinned; distance in ith bin))}\cr
     \cr
     The sequence is 1a, 2a, ..., 9a, 1b, ..., 9f.
     \cr\cr
     \item
     \code{VAR} : (\code{length(bin)-1})\eqn{\times}{ x }(\eqn{m(m-1)s^2/2})
     matrix.
     function \eqn{V} if \code{m=1}. Otherwise, the
     covariances are returned in the following ordering.
     The index for
     the lower triangle of the covariance matrix (including the
     diagonal) runs fastest, then the index of the
     conditioning species, then the index for the species the marks
     belong to. 
     That is,
     \cr
     \eqn{Cov(M_1(S_1),M_1(S_1)|d_{11}=r)}{Cov[M_1(S_1),M_1(S_1) | d_11=r]},
     \eqn{Cov(M_1(S_1),M_2(S_1)|d_{11}=r)}{Cov[M_1(S_1),M_2(S_1) | d_11=r]}, ...,
     \eqn{Cov(M_1(S_1),M_m(S_1)|d_{11}=r)}{Cov[M_1(S_1),M_m(S_1) | d_11=r]}, 
     \eqn{Cov(M_2(S_1),M_2(S_1)|d_{11}=r)}{Cov[M_2(S_1),M_2(S_1) | d_11=r]}, ...,
     \eqn{Cov(M_2(S_1),M_m(S_1)|d_{11}=r)}{Cov[M_2(S_1),M_m(S_1) | d_11=r]}, ...,
     \eqn{Cov(M_m(S_1),M_m(S_1)|d_{11}=r)}{Cov[M_m(S_1),M_m(S_1) | d_11=r]}, ...,
     \cr
     \eqn{Cov(M_1(S_1),M_1(S_1)|d_{12}=r)}{Cov[M_1(S_1),M_1(S_1) | d_12=r]}, ...,
     \eqn{Cov(M_m(S_1),M_m(S_1)|d_{12}=r)}{Cov[M_m(S_1),M_m(S_1) | d_12=r]}, ...,
     \cr
     \eqn{Cov(M_1(S_1),M_1(S_1)|d_{1s}=r)}{Cov[M_1(S_1),M_1(S_1) | d_1s=r]}, ...,
     \eqn{Cov(M_m(S_1),M_m(S_1)|d_{1s}=r)}{Cov[M_m(S_1),M_m(S_1) | d_1s=r]},
     \cr...,
     \cr
     \eqn{Cov(M_1(S_2),M_1(S_2)|d_{21}=r)}{Cov[M_1(S_2),M_1(S_2) | d_21=r]}, ...,
     \eqn{Cov(M_m(S_2),M_m(S_2)|d_{21}=r)}{Cov[M_m(S_2),M_m(S_2) | d_21=r]}, ...,
     \cr
     \eqn{Cov(M_1(S_2),M_1(S_2)|d_{2s}=r)}{Cov[M_1(S_2),M_1(S_2) | d_2s=r]}, ...,
     \eqn{Cov(M_m(S_2),M_m(S_2)|d_{2s}=r)}{Cov[M_m(S_2),M_m(S_2) | d_2s=r]},
     \cr..., 
     \cr
     \eqn{Cov(M_m(S_s),M_m(S_s)|d_{ss}=r)}{Cov[M_m(S_s),M_m(S_s) | d_ss=r]}.
     \item
     \code{VARTest} : matrix of 66 rows and the same structure of columns as
     \code{VAR}.  See \code{VAR} and \code{ETest} for details.
     \item
     \code{SQ} : (\code{length(bin)-1})\eqn{\times}{ x }(\eqn{m(m-1)s^2/2})
     matrix.  Equals sign\eqn{(V) \sqrt{|V|}}{(V) * sqrt(|V|)}.
     \item
     \code{SQTest} : Test results for SQ.  Same matrix
     dimensions as \code{VARTest}.
     \item
     \code{KMM} : Stoyan's \eqn{k_{mm}}{kmm} function\cr
     (\code{length(bin)-1})\eqn{\times}{ x }(\eqn{ms(ms+1)/2}) matrix. 
     The columns of \code{KMM} are returned in the following ordering:
     It is the lower triangle (including the diagonal) of
     expectation of the matrix \eqn{vv^\top}{v v^T}.  Here \eqn{v} is
     vector where first \eqn{m} components are the marks of species
     \eqn{S_1} the next components are the marks of \eqn{S_2} and so on.
     That is,
     \cr
     \eqn{KMM(M_1(S_1),M_1(S_1)|d_{11}=r)},...,
     \eqn{KMM(M_m(S_1),M_1(S_1)|d_{11}=r)},
     \eqn{KMM(M_1(S_2),M_1(S_1)|d_{21}=r)},...,
     \eqn{KMM(M_m(S_2),M_1(S_1)|d_{21}=r)},...,
     \eqn{KMM(M_1(S_s),M_1(S_1)|d_{s1}=r)},...,
     \eqn{KMM(M_m(S_s),M_1(S_1)|d_{s1}=r)},
     \cr...,\cr
     \eqn{KMM(M_2(S_1),M_2(S_1)|d_{11}=r)},...,
     \eqn{KMM(M_m(S_1),M_2(S_1)|d_{11}=r)},
     \eqn{KMM(M_1(S_2),M_2(S_1)|d_{21}=r)},...,
     \eqn{KMM(M_m(S_2),M_2(S_1)|d_{21}=r)},...,
     \eqn{KMM(M_1(S_s),M_2(S_1)|d_{s1}=r)},...,
     \eqn{KMM(M_m(S_s),M_2(S_1)|d_{s1}=r)},
     \cr...,\cr     
     \eqn{KMM(M_m(S_1),M_m(S_1)|d_{11}=r)},
     \eqn{KMM(M_1(S_2),M_m(S_1)|d_{21}=r)},...,
     \eqn{KMM(M_m(S_s),M_m(S_1)|d_{s1}=r)},
     \cr\cr
     \eqn{KMM(M_1(S_2),M_1(S_2)|d_{22}=r)},...,
     \eqn{KMM(M_m(S_2),M_1(S_2)|d_{22}=r)},
     \eqn{KMM(M_1(S_3),M_1(S_2)|d_{32}=r)},...,
     \eqn{KMM(M_m(S_3),M_1(S_2)|d_{32}=r)},...,
     \eqn{KMM(M_1(S_s),M_1(S_2)|d_{s2}=r)},...,
     \eqn{KMM(M_m(S_s),M_1(S_2)|d_{s2}=r)},
     \cr...,\cr
     \eqn{KMM(M_2(S_2),M_2(S_2)|d_{22}=r)},...,
     \eqn{KMM(M_m(S_2),M_2(S_2)|d_{22}=r)},
     \eqn{KMM(M_1(S_3),M_2(S_2)|d_{32}=r)},...,
     \eqn{KMM(M_m(S_3),M_2(S_2)|d_{41}=r)},...,
     \eqn{KMM(M_1(S_s),M_2(S_2)|d_{s2}=r)},...,
     \eqn{KMM(M_m(S_s),M_2(S_2)|d_{s2}=r)},
     \cr...,\cr     
     \eqn{KMM(M_m(S_2),M_m(S_2)|d_{22}=r)},
     \eqn{KMM(M_1(S_3),M_m(S_2)|d_{21}=r)},...,
     \eqn{KMM(M_m(S_s),M_m(S_2)|d_{s2}=r)},
     \cr...,\cr
     \eqn{KMM(M_1(S_s),M_1(S_s)|d_{ss}=r)},...,
     \eqn{KMM(M_m(S_s),M_1(S_s)|d_{ss}=r)},
     \eqn{KMM(M_2(S_s),M_2(S_s)|d_{ss}=r)},...,
     \eqn{KMM(M_m(S_s),M_2(S_s)|d_{ss}=r)},...,
     \eqn{KMM(M_m(S_s),M_m(S_s)|d_{ss}=r)},

     \item
     \code{GAM} : mark variogramm\cr
     (\code{length(bin)-1})\eqn{\times}{ x }(\eqn{ms(ms+1)/2}) matrix;
     see KMM for details.
     \item
     \code{Ebin} : number of values a binned value of \eqn{E} is based on.
     \item
     \code{VARbin} : number of values a binned value of \eqn{V} or
     \eqn{SQ} is based on. 
     \item
     \code{KMMbin} : number of values a binned value of
     Stoyan's \eqn{k_{mm}}{kmm} function is based on.
     \item
     \code{GAMbin} : number of values a binned value of the mark variogram
     is based on.
     \item
     \code{midbin} : centers of the bins
     \item
     \code{call}   : \code{match.call()} of \command{mpp.characteristics}
   }
  \code{summarize=FALSE} : the column structure of the above
  matrices is \code{rep} times repeated.
}
\value{
  \command{mpp.characteristics} returns
  \code{list(E, ETest, VAR, VARTest, SQ, SQTest, KMM, GAM, Ebin, 
    VARbin, KMMbin, GAMbin, midbin, call = match.call())};
  see Details. The return is \command{invisible} if \code{show=TRUE}.
}
\references{
  Schlather, M., Ribeiro, P., and Diggle, P. (2004)
  Detecting Dependence Between Marks and Locations of Marked Point Processes
  \emph{J. R. Statist. Soc., Ser. B} \bold{}, .

  Schlather, M. (2001) On the second order characteristics of marked
  point processes, \emph{Bernoulli} \bold{7}, 99-117. 
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\seealso{\command{\link{rfm.test}}, \command{\link{simulateMPP}} }

\examples{
%    library(MarkedPointProcess, lib="~/TMP")
data(BITOEK)

bin <- c(-1, seq(0, 50, 2))
normalize <- TRUE
mpp <- mpp.characteristics(bin=bin, normalize=normalize, show=TRUE,
                   coord=steigerwald$coord, diam=steigerwald$diam)
str(mpp)
}
\keyword{spatial}% at least one, from doc/KEYWORDS



