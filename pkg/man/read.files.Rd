\name{read.files}
\alias{read.files}
\title{
Read Deletion Data and Map File (if available)
}
\description{
Reads in deletion data and map file. 
}
\usage{
read.files(datafile, mapfile,   na.strings="NA", mrks.as.rows=TRUE, row.names=1, names.pres=TRUE, sep="", ...)

}
\arguments{
  \item{datafile}{
          the name of the file which the data are to be read from.
          If an  absolute path is not supplied, then the file name is
          relative to the current working directory, \code{getwd()}.
}
  \item{mapfile}{
          the name of the file containing the marker map.  
          If an  absolute path is not supplied, then the file name is
          relative to the current working directory, \code{getwd()}. If no 
          \code{mapfile} is specified, then all the marker loci are assumed
          to belong to the same chromosome.
}
  \item{sep}{
the field separator character.  Values on each line of the
          files are separated by this character.  If \code{sep = ""} (the
          default for \code{read.files}) the separator is \code{white space},
          that is one or more spaces, tabs, newlines or carriage
          returns. Both files are assumed to have the same field separator.
}
  \item{na.strings}{
a character vector of strings which are to be interpreted
          as \code{NA} values in the data file.
}
  \item{mrks.as.rows}{
a logical value indicating whether rows of the data file are the markers (\code{TRUE}) 
or the plants (\code{FALSE}).
}
  \item{row.names}{
an integer value giving the column number of the marker or plant names in the data file. 
}
  \item{names.pres}{
a logical value indicating if the data file contains the marker and plant names. If \code{FALSE}, 
\code{row.names} is ignored. 
}
  \item{...}{
other allowable \code{link{read.table}} options.
}
}
\details{
\code{read.files} is a general function for reading in deletion data and the map file (if known). 

When reading in the deletion data, 
\code{read.files} default behaviour is 
to read in a matrix of values where the rows are the markers and the columns are the plants, the 
marker and plant names are present, the first column contains the marker labels, it is a space 
seperated file, and missing genotypes are denoted by \code{NA}. 

This default behaviour can be changed via the about function arguments. 
For example, suppose the data 
were organized by rows as plants and columns as markers. Also, suppose the file is comma 
seperated and missing genotypes are denoted by 99. Then the command would be 
\code{read.files(file, mrks.as.rows=FALSE, sep=",", na.strings = "99")}.

The data can be integer or character. The integer/character value with the lowest frequency is 
assumed to denote a deletion. The values are converted into 0/1 data where 1 denotes a deletion.
 
If \code{names.pres=FALSE}, generic plant names L1, L2, \ldots and
marker names M1, M2, \ldots will be assigned to the rows and columns.  

When reading in a map file, it is assumed that it contains  three columns.
The first column has the marker names, the second the chromosome names/numbers, and the third the
cM positions of the marker loci on the chromosomes. The marker loci do not have to be in map
order. The first row must contain the column headings but these headings can be set by the 
user.  The first two columns can contain alphanumeric
values but the last column must only contain numeric values.

Missing values are not allowed and will
cause an error.

As an example of what the map file may look like:
\tabular{ccc}{
Mrk  \tab Chrm \tab  Pos  \cr
snp254  \tab 1A  \tab 23.5 \cr
snp1832 \tab 1A  \tab 104.3 \cr
snp111  \tab 2B  \tab 0.5  \cr
snp9623 \tab 1A  \tab 11.5 \cr
snp132  \tab 1A  \tab 83.5 \cr
snp018  \tab 2B  \tab 35.4 \cr
}


}
\value{
An  \code{\link{dmdata}} object is returned. 
}
\author{
Andrew W. George
}
\note{
If \code{mapfile} is not specified, then all marker loci are assumed to belong 
to the same chromsome and are mapped to the same cM position.}
\seealso{
See Also \code{\link{sim.dmdata}}, \code{\link{read.table}}, \code{\link{dmdata}}
}
\examples{
# read in deletion data
package.path = system.file("data", package="delmap")
datafn <- paste(package.path, "/deldata.txt", sep="")
mapfn <- paste(package.path, "/map.txt", sep="")

res <- read.files(datafile=datafn, mapfile = mapfn, mrks.as.rows=FALSE)

res[["x"]]


}
