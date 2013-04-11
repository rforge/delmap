\name{dmdata}
\alias{dmdata}
\alias{as.dm.matrix}
\alias{as.dm}
\alias{is.dm}
\title{Delmap.data}
\description{
\code{dmdata} is an object of  class \code{dmdata}.

\code{is.dm} tests if its argument is of object class \code{dmdata}.

\code{as.dm} attempts to turn its argument into an object of class \code{dmdata}.
}
\usage{
is.dm(dmat)

\method{as.dm}{matrix}(x)
}
\arguments{
 \item{dmat}{an object of class \code{dmdata}.
  }
 \item{x}{a \code{matrix} object.
  } 
}
\details{
\code{is.dm}  returns TRUE if x is an object of class \code{dmdata}.

\code{as.dm}  turns an object of class \code{matrix} into and object of class \code{dmdata}. 
Here, elements in the matrix are turned into 0,1 for non-deletion and deletion, respectively. Also 
generic plant and marker names are added to the rows and columns of the matrix, respectively. If the matrix 
is of the wrong form, the function returns an error message. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{ReadData}}
}


