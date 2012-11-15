\name{ReadData}
\alias{ReadData}
\title{
Read Deletion Data
}
\description{
Reads a file in table format and creates a delmap.data object from it, with 
plants corresponding to rows and genetic markers corresponding to columns 
in the file
}
\usage{
ReadData(datafile, line.names = FALSE, na.strings = "NA", marker.names = FALSE, csv = FALSE)
}
\arguments{
  \item{datafile}{
          the name of the file which the data are to be read from.
          Each row of the table appears as one line of the file.  If it
          does not contain an absolute path, the file name is
          relative to the current working directory, 'getwd()'.
}
  \item{line.names}{
a logical value indicating whether the file contains plant names as the 
first column. 
}
  \item{na.strings}{
a character vector of strings which are to be interpreted
          as 'NA' values. 
}
  \item{marker.names}{
 a logical value indicating whether the file contains the
          names of the marker loci as its first line.  
}
  \item{csv}{
a logical value indicating if the file is comma seperated. 
}
}
\details{
Deletion data are read into R with this function. The data can be integer or character where any integer/character value is valid. The integer/character
value with the lowest frequency is assumed to denote a deletion. The values are converted into 0/1 data where 1 denotes a deletion. 
If row/column names are not specified in this file, plant names L1, L2, \ldots and 
marker names M1, M2, \ldots will be assigned to the rows and columns, respectively.  
}
\value{
An object of class delmap.data. It contains the data in matrix form that has 0/1 elements. It also contains the plant 
names and marker loci names. 
}
\author{
Andrew W. George
}
\seealso{
See Also \code{\link{SimDeletions}}
}
\examples{
\dontrun{
# read in deletion data 
deldat <- ReadData(datafile="./deldata.txt", marker.names=TRUE, line.names=TRUE)
}
}
