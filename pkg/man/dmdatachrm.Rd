\name{dmdatachrm}
\alias{dmdatachrm}
\alias{as.dmdatachrm.matrix}
\alias{as.dmdatachrm}
\alias{is.dmdatachrm}
\title{dmdatachrm}
\description{
\code{dmdatachrm} is an object of  class \code{dmdatachrm} which is 
escentially a \code{matrix} 
object with \code{attr(dmdata, "dimnames")[1]} being a character vector of 
the line names, 
\code{attr(dmdata, "dimnames")[1]} being a character vector 
of the marker names, and 
\code{attr(dmdata, "imputed")} being an additional attribute which is 
a two columned matrix giving the row,column 
location of any elements in \code{dmdatachrm} that have had 
their values imputed.
Unlike the object class \code{\link{dmdata}}, there is no \code{nchrm} 
attribute because all marker data in this class is assumed to 
belong to the same 
chromosome. 


\code{is.dmdatachrm} tests if its argument is of object class \code{dmdatachrm}.

\code{as.dmdatachrm} attempts to turn its argument into an 
object of class \code{dmdatachrm}.
}
\usage{
is.dmdatachrm(dmat)

\method{as.dmdatachrm}{matrix}(x)
}
\arguments{
 \item{dmat}{an object of class \code{dmdatachrm}.
  }
 \item{x}{a \code{matrix} object.
  } 
}
\details{
\code{is.dmdatachrm}  returns TRUE if x is an object of class \code{dmdatachrm}.

\code{as.dmdatachrm}  turns an object of class \code{matrix} into and object 
of class \code{dmdatachrm}. It is assumed that 
\itemize{
\item the rows of the matrix are the plants and the columns of the matrix 
are the marker loci
\item the marker loci belong to the same chromosome. 
}
The  elements of the matrix are turned into 0,1 for 
non-deletion and deletion, respectively. Also 
generic plant and marker names are added to the rows and columns 
of the matrix, respectively. If the matrix 
is of the wrong form, the function returns an error message. 
}
\author{
Andrew W. George
}
\seealso{
\code{\link{read.files}}
}


