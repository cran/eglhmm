\encoding{UTF-8}
\name{weissData}
\alias{weissData}
\alias{Bovine}
\alias{Cryptosporidiosis}
\alias{Downloads}
\alias{EricssonB_Jul2}
\alias{FattyLiver}
\alias{FattyLiver2}
\alias{goldparticle380}
\alias{Hanta}
\alias{InfantEEGsleepstates}
\alias{IPs}
\alias{LegionnairesDisease}
\alias{OffshoreRigcountsAlaska}
\alias{PriceStability}
\alias{Strikes}
\alias{WoodPeweeSong}
\docType{data}
\title{
   Data from \dQuote{An Introduction to Discrete-Valued Time Series}
}
\description{
   Data sets from the book
   \dQuote{An Introduction to Discrete-Valued Time Series} by
   Christian H. \enc{Weiß}{Weiss}.
%% To get the ß symbol use <compose key> s s.

   The data sets are named \code{Bovine}, \code{Cryptosporidiosis},
   \code{Downloads}, \code{EricssonB_Jul2}, \code{FattyLiver},
   \code{FattyLiver2}, \code{goldparticle380},
   \code{Hanta}, \code{InfantEEGsleepstates}, \code{IPs},
   \code{LegionnairesDisease}, \code{OffshoreRigcountsAlaska},
   \code{PriceStability}, \code{Strikes} and \code{WoodPeweeSong}.
}
\format{
  Each data set is a data frame with a single column named \code{"y"}.
  \itemize{
  \item \code{Bovine}
  There are 8419 rows.  The column \code{"y"} is a factor,
  with levels \code{"a","c","g","t"}, the DNA \dQuote{bases}.
  It constitutes the DNA sequence of the bovine leukemia virus.
  \item \code{Cryptosporidiosis}
  There are 365 rows.  The column \code{"y"} is a numeric (integer)
  vector.  It consists of weekly counts of new infections, in Germany
  in the years 2002 to 2008.  The counts vary between 2 and 78.
  \item \code{Downloads}
  There are 267 rows.  The column \code{"y"} is a numeric
  (integer) vector.  It consists of the daily number of downloads
  of a TEX editor for the period from June 2006 to February 2007.
  These counts vary between 0 and 14.
  \item \code{EricssonB_Jul2}
  There are 460 rows.  The column \code{"y"} is a numeric (integer)
  vector.  It consists of the number of transactions per minute,
  of the Ericsson B stock, between 9:35 and 17:14 on 2 July, 2002.
  The counts vary between 0 and 37.
  \item \code{FattyLiver}
  There are 928 rows.  The column \code{"y"} is a numeric (binary)
  vector.  The value 1 indicates that \dQuote{the considered
  diagnosis cannot be excluded for the current patient; that is,
  suitable countermeasures are required}, and the value 0 indicates
  that this is not so.  The values refer to different patients,
  examined sequentially over time.
  \item \code{FattyLiver2}
  There are 449 rows.  The column \code{"y"} is a numeric (binary)
  vector as for \code{FattyLiver}.  (Different examiner, different
  sequence of patients.)
  \item \code{goldparticle380}
  There are 380 rows.  The column \code{"y"} is a numeric (integer)
  vector of counts of gold particles measured in a fixed volume
  element of a colloidal solution over time. The count values
  vary because of the Brownian motion of the particles.  They vary
  between 0 and 7.
  \item \code{Hanta}
  There are 52 rows.  The column \code{"y"} is a numeric (integer)
  vector consisting of the weekly number of territorial units (out of \eqn{n = 38}
  territorial units with at least one new case of a hantavirus infections,
  in the year 2011.  The numbers vary between 0 and 11. 
  \item \code{InfantEEGsleepstates}
  There are 107 rows.  The column \code{"y"} is a factor with
  levels \code{qt, qh, tr, al, ah, aw}.  The level \code{"aw"}
  does not actually appear.
  \item \code{IPs}
  There are 241 rows.  The column \code{"y"} is a numeric (integer)
  vector of the counts of different IP addresses registered at a
  web server within periods of length two minutes, \dQuote{assumed}
  to have been observed between 10:00 a.m. and 6:00 p.m. on 29
  November 2005.  The counts vary between 0 and 8.
  \item \code{LegionnairesDisease}
  There are 365 rows.  The column \code{"y"} is a numeric (integer)
  vector of weekly counts of new infections in Germany, in the
  years 2002 to 2008.  The counts vary between 0 and 26.
  \item \code{OffshoreRigcountsAlaska}
  There are 417 rows.  The column \code{"y"} is a numeric (integer)
  vector of weekly counts of active rotary drilling rigs in Alaska
  for the period 1990 to 1997. The counts vary between 0 and 6.
  \item \code{PriceStability}
  There are 152 rows.  The column \code{"y"} is a numeric (integer)
  vector of monthly counts of countries (out of a group of 17
  countries) that showed stable prices (that is, an inflation rate
  below 2\%), in the period from January 2000 to December 2006.
  The counts vary between 0 and 17.
  \item \code{Strikes}
  There are 108 rows.  The column \code{"y"} is a numeric (integer)
  vector of the monthly counts of work stoppages (strikes and
  lock-outs) of 1000 or more workers in the period 1994 to 2002.
  The counts vary between 0 and 14.
  \item \code{WoodPeweeSong}
  There are 1327 rows.  The column \code{"y"} is a factor with
  levels \code{"1", "2", "3"} corresponding to the three different
  \dQuote{phrases} of wood wewee song.  The time series comprises
  a sequence of observations of the \dQuote{morning twilight}
  song of the wood pewee.
  }
}

\details{
   For detailed information about each of these data sets, see the
   book cited in the \bold{References}.

   Note that the data sets \code{Cryptosporidiosis}
   and \code{LegionnairesDisease} are actually
   called \cr \code{Cryptosporidiosis_02-08} and
   \code{LegionnairesDisease_02-08} in the given reference.
   The \cr \dQuote{suffixes} were removed since the minus sign causes
   problems in a variable name in \code{R}.
}

\source{
These data sets were kindly provided by Prof. Christian
H. \enc{Weiß}{Weiss}.  The package author is also pleased
to acknowledge the kind permission granted by Prof. Kurt
\enc{Brännäs}{Brannas} (Professor Emeritus of Economics at
\enc{Umeå}{Umea} University) to include the Ericsson time series
data set (\code{EricssonB_Jul2}).
}
%% Email: weissc@hsu-hh.de
\references{Christian H. \enc{Weiß}{Weiss} (2018).  \emph{An
Introduction to Discrete-Valued Time Series}.  Chichester: John
Wiley & Sons.
}
\examples{
\dontrun{
fit1 <- hmm(WoodPeweeSong,K=2,verbose=TRUE)
# EM converges in 6 steps --- suspicious.
set.seed(321)
fit2 <- hmm(WoodPeweeSong,K=2,verbose=TRUE,rand.start=list(tpm=TRUE,Rho=TRUE))
# 52 steps --- note the huge difference between fit1$log.like and fit2$log.like!
set.seed(321)
fit3 <- hmm(WoodPeweeSong,K=2,verbose=TRUE,method="bf",
            rand.start=list(tpm=TRUE,Rho=TRUE))
# log likelihood essentially the same as for fit2
}
}
\keyword{datasets}
