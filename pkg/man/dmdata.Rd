\name{dmdata}
\alias{dmdata}
\alias{as.dmdata.matrix}
\alias{as.dmdata.list}
\alias{as.dmdata.dmdatachrm}
\alias{as.dmdata}
\alias{is.dmdata}
\title{dmdata}
\description{
\code{dmdata} is an object of  class \code{dmdata} which is 
escentially a \code{list} 
object. Each element of the list is a \code{\link{dmdatachrm}}
object. The number of elements in the list corresponds to the number 
of chromosomes with marker data. The \code{dmdata} object has 
an added attribute \code{attr(dmdata, "nchrm")} for the number of 
chromosomes with data.   

\code{is.dmdata} tests if its argument is of object class \code{dmdata}.

\code{as.dmdata} attempts to turn its argument into an object of class \code{dmdata}.
}
\usage{
is.dmdata(x)

\method{as.dmdata}{dmdatachrm}(x)

\method{as.dmdata}{matrix}(x)

\method{as.dmdata}{list}(x)
}

\arguments{
 \item{x}{an object of class \code{dmdata} when used with \code{is.dmdata}.
  }
 \item{x}{an object of type \code{dmdatachrm}, \code{matrix}, or \code{list}
  } 
}


\details{
\code{is.dmdata}  returns TRUE if x is an object of class \code{dmdata}.

\code{as.dmdata}  turns an object of class \code{matrix} 
or \code{\link{dmdatachrm}} into and object 
of class \code{dmdata}. If the object is a \code{matrix}, then 
it is assumed that 
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
\code{\link{as.dmdatachrm}}, \code{\link{is.dmdatachrm}}
}


