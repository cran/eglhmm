\name{plot.eglhmm}
\alias{plot.eglhmm}
\title{
    Plot the results of an eglhmm fit.
}
\description{
   For each specified model cell plot an array, with one panel
   for each state, of the probability mass or density functions
   corresponding to the given cell and state.  The plots are
   produced with \code{type="h"} for probability mass functions,
   and with \code{type="l"} for probability density functions.
}
\usage{
   \method{plot}{eglhmm}(x, ..., wcells = NULL, col = "red",
             nrnc = NULL, ntop = NULL, xlab = NULL, ylab = NULL,
             xlim = NULL, ylim = NULL, main = NULL, cex.main = 1.5)
}
\arguments{
  \item{x}{
  An object of class \code{"eglhmm"} as returned by the function
  \code{\link{eglhmm}()}.
}
  \item{\dots}{
  Not used.
}
  \item{wcells}{
  Character vector specifying the cells of the model to be plotted.
  Defaults to all cells (i.e. the levels of \code{x$data$cf}).
}
  \item{col}{
  The colour for the (vertical) lines of the plots.
}
  \item{nrnc}{
  An integer vector of length two specifying the dimenions of the
  array of plots that is produced.  The first entry is the number
  of rows, the second the number of columns.  The product of the
  entries must be greater than or equal to \code{K}, the number of
  states in the model.
}
  \item{ntop}{
  The largest \code{x}-value to be used in plots of the Poisson
  distribution.  Defaults
  to the maximum of the upper \eqn{10^{-7}}{1e-7} quantile of
  all of the Poisson distributions that are to be plotted.
  Ignored unless \code{x$distr} is \code{"Poisson"}.
}
  \item{xlab}{
  An optional label for the \eqn{x}-axes of the panels in the array
  of plots.  Defaults to \code{"x"}.
}
  \item{ylab}{
  An optional label for the \eqn{y}-axes of the panels in the array
  of plots.  Defaults to \code{"probability"} for probability mass
  functions and to \code{"probability density"} for probability
  density functions.
}
  \item{xlim}{
  An optional vector of length two, specifying the \code{x}-limits
  for the plots.  Defaults to \code{c(0,ntop)} if \code{x$distr} is
  \code{"Poisson"} and to \code{c(0,x$size)} if \code{x$distr} is
  \code{"Binomial"}.  There is no default if \code{x$distr} is
  \code{"Gaussian"}.
}
  \item{ylim}{
  An optional vector of length two, specifying the \code{y}-limits
  for the plots.  Defaults to \code{c(0,M)}  where \code{M} is
  the maximum of all of the probabilities or probability density
  values that are to be plotted.
}
  \item{main}{
  Optional character vector specifying overall titles for each
  array panel of plots.  Defaults to the names of the model cells.
  If the length of \code{main} is less than the number (\code{nwc})
  of cells to be plotted, then main is replicated to have length
  \code{nwc}.  If the length of main is greater than \code{nwc}
  then entries with index greater than \code{nwc} are ignored.
  If you wish there to be no overall titles for the arrays of plots,
  specify \code{main=""}.
}
  \item{cex.main}{
  Expansion factor for the text in the main title (determining
  the size of the text).  Ignored if \code{main} is set equal
  to the empty string.
}
}
\details{
  If plotting is interactive, then the arrays of plots are displayed
  one at a time, and (except for the last of the plots) the user is
  prompted with the string \code{"Go?"} after each array is plotted.
  Press \code{<return>} to see the next plot.
}
\value{
  None.
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}

\seealso{
  \code{\link{eglhmm}()}
}
\examples{
    loc4 <- c("LngRf","BondiE","BondiOff","MlbrOff")
    SCC4 <- SydColCount[SydColCount$locn \%in\% loc4,]
    SCC4$locn <- factor(SCC4$locn) # Get rid of unused levels.
    rownames(SCC4) <- 1:nrow(SCC4)
    fit <- eglhmm(y~locn+depth,data=SCC4,cells=c("locn","depth"),
                 K=2,distr="P",verb=TRUE)
    plot(fit)
    allcells <- levels(fit$data$cf)
    wcells   <- allcells[grep("\\\\.60",allcells)]
    plot(fit,wcells=wcells,main=c("Longreef","Bondi East","Bondi Offshore",
                                  "Malabar Offshore"),ntop=12)
}
\keyword{ hplot }
