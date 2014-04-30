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
object. The elements of \code{dmdata} are the deletion data for a chromosome and the  marker map. 
The deletion data for a chromosome is a  \code{\link{dmdatachrm}} object.  
The \code{dmdata} object has 
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
, \code{\link{dmdatachrm}}, or \code{list} into and object 
of class \code{dmdata}. If the object is a \code{matrix}, then 
it is assumed that 
\itemize{
\item the rows of the matrix are the plants and the columns of the matrix 
are the marker loci
\item the marker loci belong to the same chromosome. 
}
Here, the  elements of the matrix are turned into 0,1 for 
non-deletion and deletion, respectively. Also 
generic plant and marker names are added to the rows and columns 
of the matrix, respectively.  A generic map is also generated where the marker loci are assumed to be from a single 
chromosome. 
If the matrix 
is of the wrong form, the function returns an error message. 

If the object is of class \code{\link{dmdatachrm}}, then it is assumed that
all the marker loci are located on the same chromosome at the same cM position. 


If it is a \code{list} object, then it is assumed that 
\itemize{
\item the list contains a \code{map} element
\item the non-map elements are objects of class \code{\link{dmdatachrm}} 
}




}
\author{
Andrew W. George
}
\seealso{
\code{\link{as.dmdatachrm}}, \code{\link{is.dmdatachrm}}
}


