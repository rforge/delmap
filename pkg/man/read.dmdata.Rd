\name{read.dmdata}
\alias{read.dmdata}
\title{
Read Deletion Data
}
\description{
Reads in deletion data. 
}
\usage{
read.dmdata(file, sep="", na.strings="NA", mrks.as.rows=TRUE, row.names=1, names.pres=TRUE, ...)
}
\arguments{
  \item{file}{
          the name of the file which the data are to be read from.
          If an  absolute path is not supplied, then the file name is
          relative to the current working directory, \code{getwd()}.
}
  \item{sep}{
the field separator character.  Values on each line of the
          file are separated by this character.  If \code{sep = ""} (the
          default for \code{read.dmdata}) the separator is \code{white space},
          that is one or more spaces, tabs, newlines or carriage
          returns.
}
  \item{na.strings}{
a character vector of strings which are to be interpreted
          as \code{NA} values.
}
  \item{mrks.as.rows}{
a logical value indicating whether rows of the file are the markers (\code{TRUE}) 
or the plants (\code{FALSE}).
}
  \item{row.names}{
an integer value giving the column number of the marker or plant names. 
}
  \item{names.pres}{
a logical value indicating if the file contains the marker and plant names. If \code{FALSE}, 
\code{row.names} is ignored. 
}
  \item{...}{
other allowable \code{link{read.table}} options.
}
}
\details{
\code{read.dmdata} is a general function for reading in deletion data. Its default behaviour is 
to read in a matrix of values where the rows are the markers and the columns are the plants, the 
marker and plant names are present, the first column contains the marker labels, it is a space 
seperated file, and missing genotypes are denoted by \code{NA}. 

This default behaviour can be changed via the about function arguments. 
For example, suppose the data 
were organized by rows as plants and columns as markers. Also, suppose the file is comma 
seperated and missing genotypes are denoted by 99. Then the command would be 
\code{read.dmdata(file, mrks.as.rows=FALSE, sep=",", na.strings = "99")}.

The data can be integer or character. The integer/character value with the lowest frequency is 
assumed to denote a deletion. The values are converted into 0/1 data where 1 denotes a deletion.
 
If \code{names.pres=FALSE}, generic plant names L1, L2, \ldots and
marker names M1, M2, \ldots will be assigned to the rows and columns.  

Any valid \code{\link{read.table}} can be used. 

}
\value{
An object of class \code{dmdata}. It contains the data in matrix form that has 0/1 as its 
elements. Missing genotypes are (re)coded as \code{NA}.
Irregardless of the value of \code{mrks.as.rows}, an object is returned where the 
rows are the markers and the columns are the plants. 
}
\author{
Andrew W. George
}
\seealso{
See Also \code{\link{sim.dmdata}}, \code{\link{read.table}}, \code{\link{dmdata}}
}
\examples{
# read in deletion data
package.path = system.file("data", package="delmap")
filename <- paste(package.path, "/deldata.txt", sep="")

deldat <- read.dmdata(file=filename, mrks.as.rows=FALSE)

deldat


}
