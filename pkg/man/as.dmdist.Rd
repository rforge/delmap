\name{dmdist}
\alias{dmdist}
\alias{as.dmdist}
\alias{is.dmdist}
\alias{as.dmdist.list}
\title{Distance}
\description{
\code{dmdist} is a multivariate distance object of  class \code{dmdata}.

\code{is.dmdist} tests if its argument is of object class \code{dmdata}.

\code{as.dmdist} attempts to turn its argument into an object of class \code{dmdata}.
}
\usage{
is.dmdist(y)

\method{as.dmdist}{list}(x)
}
\arguments{
 \item{y}{an object of class \code{dmdist}.
  }
 \item{x}{a \code{list} object.
  } 
}
\details{
\code{is.dmdist}  returns TRUE if x is an object of class \code{dmdist}.

\code{as.dmdist}  turns an object of class \code{list} into and object of class \code{dmdist}. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{RearrangeDist}}
}


