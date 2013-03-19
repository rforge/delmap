\name{delmap.data}
\alias{delmap.data}
\alias{as.delmap.matrix}
\alias{as.delmap}
\alias{is.delmap}
\title{Delmap.data}
\description{
\code{delmap.data} is an object of  class \code{delmap.data}.

\code{is.delmap} tests if its argument is of object class \code{delmap.data}.

\code{as.delmap} attempts to turn its argument into an object of class \code{delmap.data}.
}
\usage{
is.delmap(dmat)

\method{as.delmap}{matrix}(x)
}
\arguments{
 \item{dmat}{an object of class \code{delmap.data}.
  }
 \item{x}{a \code{matrix} object.
  } 
}
\details{
\code{is.delmap}  returns TRUE if x is an object of class \code{delmap.data}.

\code{as.delmap}  turns an object of class \code{matrix} into and object of class \code{delmap.data}. 
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


