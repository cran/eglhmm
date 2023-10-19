\name{SydColDat}
\alias{SydColCount}
\alias{SydColDisc}
\docType{data}
\title{Sydney coliform bacteria data}
\description{
   Transformed counts of faecal coliform bacteria in sea water
   at seven locations: Longreef, Bondi East, Port Hacking ``50'',
   and Port Hacking ``100'' (controls) and Bondi Offshore, Malabar
   Offshore and North Head Offshore (outfalls).  At each location
   measurements were made at four depths: 0, 20, 40, and 60 meters.
}
\usage{
   SydColCount
   SydColDisc
}
\format{
  Data frames with 5432 observations on the following 6 variables.
  \describe{
    \item{\code{y}}{Transformed measures of the number of faecal
    coliform count bacteria in a sea-water sample of some specified
    volume.  The original measures were obtained by a repeated
    dilution process.

    For \code{SydColCount} the transformation used was essentially
    a square root transformation, resulting values greater than 150
    being set to \code{NA}.  The results are putatively compatible
    with a Poisson model for the emission probabilities.

    For \code{SydColDisc} the data were discretised
    using the \code{cut()} function with breaks given
    by \code{c(0,1,5,25,200,Inf)} and labels equal to
    \code{c("lo","mlo","m","mhi","hi")}.}

    Note that in the \code{SydColDisc} data there are 180 fewer
    missing values (\code{NA}s) in the \code{y} column than in
    the \code{SydColCount} data.  This is because in forming
    the \code{SydColCount} data (transforming the original data
    to a putative Poisson distribution) values that were greater
    than 150 were set equal to \code{NA}, and there were 180 such
    values.

    \item{\code{locn}}{a factor with levels \dQuote{LngRf}
    (Longreef), \dQuote{BondiE} (Bondi East), \dQuote{PH50}
    (Port Hacking 50), \dQuote{PH100}  (Port Hacking 100),
    \dQuote{BondiOff} (Bondi Offshore), \dQuote{MlbrOff} (Malabar
    Offshore) and \dQuote{NthHdOff} (North Head Offshore)}

    \item{\code{depth}}{a factor with levels \dQuote{0} (0 metres),
    \dQuote{20} (20 metres), \dQuote{40} (40 metres) and \dQuote{60}
    (60 metres).}

    \item{\code{ma.com}}{A factor with levels \code{no} and \code{yes},
    indicating whether the Malabar sewage outfall had been commissioned.}

    \item{\code{nh.com}}{A factor with levels \code{no} and \code{yes},
    indicating whether the North Head sewage outfall had been commissioned.}

    \item{\code{bo.com}}{A factor with levels \code{no} and \code{yes},
    i.ndicating whether the Bondi Offshore sewage outfall had been commissioned.}
  }
}
\section{Oh for fuck's sake!!!}{
Blah.
}

\details{
   The observations corresponding to each location-depth combination
   constitute a time series.  The sampling interval is ostensibly
   1 week; distinct time series are ostensibly synchronous.
   The measurements were made over a 194 week period.  See Turner
   et al. (1998) for more detail.
}
\source{
  Geoff Coade, of the New South Wales Environment Protection
  Authority (Australia) 
}
\references{
  T. Rolf Turner, Murray A. Cameron, and Peter J. Thomson.  Hidden
  Markov chains in generalized linear models.  Canadian J. Statist.,
  vol. 26, pp. 107 -- 125, 1998.
  
  Rolf Turner.  Direct maximization of the likelihood of a hidden
  Markov model. \emph{Computational Statistics and Data Analysis}
  \bold{52}, pp. 4147 -- 4160, 2008, doi:10.1016/j.csda.2008.01.029.

}
\examples{
# Select out a subset of four locations:
loc4 <- c("LngRf","BondiE","BondiOff","MlbrOff")
SCC4 <- SydColCount[SydColCount$locn \%in\% loc4,] 
SCC4$locn <- factor(SCC4$locn) # Get rid of unused levels.
rownames(SCC4) <- 1:nrow(SCC4)
}
\keyword{datasets}
