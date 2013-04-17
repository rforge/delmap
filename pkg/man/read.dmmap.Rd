\name{read.dmmap}
\alias{read.dmmap}
\title{
Read Map File
}
\description{
Reads in marker map
}
\usage{
read.dmmap(file, sep="",   ...)
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
          default for \code{read.dmmap}) the separator is \code{white space},
          that is one or more spaces, tabs, newlines or carriage
          returns.
}
  \item{...}{
other allowable \code{link{read.table}} options.
}
}
\details{
The map file contains  three columns. 
The first column has the marker names, the second the chromosome names/numbers, and the third the 
cM positions of the marker loci on the chromosomes. The marker loci do not have to be in map 
order. The first row must contain the column headings. The first two columns can contain alphanumeric 
values but the last column must only contain numeric values. 

The default format is for the file to be space delimated. Missing values are not allowed and will 
cause an error. 

As an example of what the map file may look like:
\tabular{ccc}{
Mrk  \tab Chrm \tab  Pos  \cr
snp254  \tab 1  \tab 23.5 \cr
snp1832 \tab 1  \tab 104.3 \cr
snp111  \tab 2  \tab 0.5  \cr
snp9623 \tab 1  \tab 11.5 \cr
snp132  \tab 1  \tab 83.5 \cr
snp018  \tab 2  \tab 35.4 \cr
}

A comma seperated file is read in by setting \code{sep=","}. 
}
\value{
A list is returned where its length is equal to the number of chromosomes. Each element of the list
is a data frame   


data frame with three columns is returned where the first column contains the marker names, 
the second column contains the chromosome names, and the third column contains the position (in cM) 
of the marker loci.  Rows are placed in  map order. 
}
\author{
Andrew W. George
}
\seealso{
See Also \code{\link{SimDeletions}}, \code{\link{read.table}}, \code{\link{dmdata}}
}
\examples{
# read in deletion data
package.path = system.file("data", package="delmap")
filename <- paste(package.path, "/map.txt", sep="")

delmap <- read.dmmap(file=filename)

delmap


}
