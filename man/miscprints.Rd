\name{miscprints}
\alias{print.RhoExpForm}
\alias{print.RhoProbForm}
\alias{print.kitty}
\title{
    Specialised print methods.
}
\description{
   Print objects of class \code{"RhoExpForm"}, \code{"RhoProbForm"}
   and \code{"kitty"}, appropriately.
}
\usage{
\method{print}{RhoExpForm}(x, ...)
\method{print}{RhoProbForm}(x, ...)
\method{print}{kitty}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{"RhoExpForm"},
  \code{"RhoProbForm"} or \code{"kitty"} respectively.
}
  \item{\dots}{
  Not used.  Present for compatibility with the generic
  \code{print()} function.
}
}
\details{
  The methods \code{print.RhoExpForm()} and
  \code{print.RhoProbForm()} are present essentially for debugging
  purposes only.  The method \code{print.kitty()} is present to
  improve the appearance of printed output from \code{eglhmm}
  when there is a \code{"message"} component of this output.
  None of these methods would normally be called by users.
}
\value{
  None.
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}
\keyword{ print }
