\name{dmdata}
\alias{dmdata}
\alias{as.dmdata.matrix}
\alias{as.dmdata}
\alias{is.dmdata}
\title{Delmap.data}
\description{
\code{dmdata} is an object of  class \code{dmdata} which is escentially a \code{matrix} 
object with \code{attr(dmdata, "dimnames")[1]} being a character vector of the line names, 
\code{attr(dmdata, "dimnames")[1]} being a character vector of the marker names, and 
\code{attr(dmdata, "imputed")} being an additional attribute which is 
a two columned matrix giving the row,column 
location of any elements in \code{dmdata} that have had their values imputed.

\code{is.dmdata} tests if its argument is of object class \code{dmdata}.

\code{as.dmdata} attempts to turn its argument into an object of class \code{dmdata}.
}
\usage{
is.dmdata(dmat)

\method{as.dmdata}{matrix}(x)
}
\arguments{
 \item{dmat}{an object of class \code{dmdata}.
  }
 \item{x}{a \code{matrix} object.
  } 
}
\details{
\code{is.dmdata}  returns TRUE if x is an object of class \code{dmdata}.

\code{as.dmdata}  turns an object of class \code{matrix} into and object of class \code{dmdata}. 
Here, elements in the matrix are turned into 0,1 for non-deletion and deletion, respectively. Also 
generic plant and marker names are added to the rows and columns of the matrix, respectively. If the matrix 
is of the wrong form, the function returns an error message. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{read.dmdata}}
}


