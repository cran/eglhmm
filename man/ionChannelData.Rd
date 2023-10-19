\name{ionChannelData}
\alias{ionChannelData}
\alias{ic25kHz_12_sgmnt1}
\alias{ic25kHz_13_sgmnt2}
\alias{ic25kHz_14_sgmnt2}
\alias{ic25kHz_15_sgmnt2}
\alias{ic50kHz_06_sgmnt2}
\alias{ic50kHz_08_sgmnt2}
\alias{ic50kHz_09_sgmnt1}
\alias{ic50kHz_10_sgmnt1}
\docType{data}
\title{
    Ion channel data
}
\description{
Time series of observations, made by means of patch clamps, of
current in picoamps, across cell membranes.
}
\section{Note}{
These data are \bold{not} immediately available in the \code{eglhmm}
package.  Their presence would cause the size of the \code{data}
directory to exceed 4.5 Mb., which is unacceptably large.
Consequently these data sets have been placed in a separate
\dQuote{data only} package called \code{ionChannelData}, which
is available from \code{github}.  This package may be obtained by
executing the command:
\preformatted{
    install.packages("ionChannelData",repos="https://rolfturner.r-universe.dev")
}
After having installed the \code{ionChannelData} package, you may
load it via \code{library(ionChannelData)} and then access the data
sets in the usual way, e.g. \code{X <- ic25kHz_12_sgmnt1}.

Alternatively (after having installed the \code{ionChannelData}
package) you may use the \code{::} syntax to access a single data
set, e.g. \code{X <- ionChannelData::ic25kHz_12_sgmnt1}.

You can access the documentation via,
e.g. \code{?ionChannelData::ionChannelData}.
}

\keyword{datasets}
