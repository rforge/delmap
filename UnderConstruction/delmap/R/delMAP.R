## Copyright (C) 2012 Andrew W. George.
## Comment: this R source code is simulate and analyses deletion type data
## Methods:  treat deletion mapping as TSP
##           solve TSP with branch and cut (as implemented in condore)
##           impute missing data with simulated annealing
## Standard file inputs:  
##         number of lines, 
##         number of markers
## Standard file outputs:  
##         best ordering, 
##         plot of shortest distance, 
##         plot of all orderings
## To Do
##   1. code to bin markers with cor=1 and no missing genotypes. 



  
####--------------------------------##
##   INTERNAL FUNCTION DECLARATION  ##
####--------------------------------##
  write.it <- function(x)
  {
    ## Internal function
    ## Purpose: to write out the data

      if( nrow(x) < 5){
         if( ncol(x) < 8)
            write.table(format(x, justify="right"), row.names=F, col.names=F, quote=F)
           else
            write.table(format(x[,1:8], justify="right"), row.names=F, col.names=F, quote=F)
       } else {
         if( ncol(x) < 8)
            write.table(format(x[1:5,], justify="right"), row.names=F, col.names=F, quote=F)
           else
            write.table(format(x[1:5,1:8], justify="right"), row.names=F, col.names=F, quote=F)
       } ## end outer if

   } ## end function write.it

  write.mrknames <- function(x)
  {
    ## Internal function
    ## Purpose: to write out marker names
    cat("\n")
    cat(" Marker names (in column data order) are: \n")
    if(length(rownames(x))  < 5)
        cat(" ", rownames(x),      "\n")
      else
        cat(" ", rownames(x)[1:5], "... \n")

  } ## end function write.mrknames


  write.linenames <- function(x)
  {
    ## Internal function
    ## Purpose:  write out line names 

    cat("\n")
    cat(" Line names (in row data order) are: \n")
    if(length(colnames(x)) < 5)
       cat(" ", colnames(x),  "\n\n")
     else
       cat(" ", colnames(x)[1:5], "... \n\n")

     cat(" The proportion of genotypes missing is", round(sum(is.na(x))/(nrow(x)*ncol(x)),2), "\n")
  }  ## end function write.linenames


 read.data <- function(file, header, sep, na.strings, row.names)
  {
     # internal function
     # purpose: to read data file and create data frame
     dat <- read.table(file=file, header=header, sep= sep, na.strings=na.strings,
             row.names=row.names)
     dat <- as.matrix(dat)  ## converts data frame into matrix
    if(length(table(dat)) > 2)
    stop(" Error:  more than two different genotypes have been detected. \n   Deletion data must contain presence/absence genotypes. ")
     return(dat)
  }  ## end function read.dat


  checkinputs.ReadData <- function(df, ln, mn)
  {
    ## Internal function
    ## Purpose:  checks for presence of data file, if line.names and marker.names 
    ##           are logical variables. 
 if(!file.exists(df))
   stop(" Input file does not exist.")
 if(!is.logical(ln))
   stop("line.names must be TRUE/FALSE only. If TRUE, the first column contains the names of the lines.")
 if(!is.logical(mn))
   stop("marker.names must be TRUE/FALSE only. If TRUE, then first row of data file contains marker names")
  } ## end function check.inputs.ReadData

  create.structure.SimDeletion <- function(nlines, nmarkers, labs)
  {
   ## Internal function
   ## Purpose:     To create data structure for SimDeletions
 
   ## create 2-D data structure 
   deldat  <- matrix(data=0, nrow=nlines, ncol=nmarkers)
   if(!is.null(labs)) 
    { if(length(labs) != nmarkers)  stop("Error: number of marker labels not of correct length")
    }  else {
     ## assign generated marker labels
     labs <- paste( rep("M", nmarkers), 1:nmarkers, sep="")
    }  ## end if !is.null
    colnames(deldat) <- labs

    labs <- paste( "L", 1:nrow(deldat), sep="")
    rownames(deldat) <- labs
    return(deldat)
   } ## end function create.structure.SimDeletion



  step1.simulation <- function(numlines, plines)
  {
    ## Internal function
    ## Purpose: step 1 of simulation function
     rnumlines <- rbinom(1, numlines, prob=plines) ## rnd number of lines 
                                                 ## to contain deletions
     indx <- NULL
     if(rnumlines > 0) chosen.lines <- sample(1:numlines, rnumlines, replace=F)  ## index of lines that
                                                 ## will contain deletions
     list(rnumlines=rnumlines, chosen.lines=chosen.lines)
   } ## end function step1.simulation



  step23.simulation <- function(rnumlines, lambda, numlines, nummarkers, chosen.lines, deldat)
  {
   ## Internal function
   ## Purpose: step2 and step3 of simulation

   ## 2 & 3
    ## create vector of deletion sizes from a poisson with averge number of del Enumdel
    delsize <- rpois(rnumlines, lambda)
    delsize[delsize>=numlines] <- rnumlines


    for(ii  in 1:length(delsize))
    {
       start.rindx <- sample( 1:(1+nummarkers- delsize[ii]), 1) 
       
       deldat[ chosen.lines[ii],  start.rindx:(start.rindx+delsize[ii] - 1)  ] <- 1

    }  ## end for ii  
    deldat
   }

  ## create new data matrix with missing gentotypes if p.missing is not null
  add.missing.deldat <- function(p.missing, deldat)
  {


    if(!is.null(p.missing))
    {
      mis.deldat <- deldat
      mat.indx <- which(deldat==0 | deldat==1, arr.ind=T)
      mis.indx <- mat.indx[ sample(1:nrow(mat.indx),
                            round(nrow(deldat)*ncol(deldat)*p.missing),replace=T),]
      mis.deldat[mis.indx] <- NA
      class(mis.deldat) <- "delmap.data"
    }

    class(deldat) <- "delmap.data"
    ifelse( is.null(p.missing),
       res <- list(nomissing=deldat),
       res <- list(nomissing=deldat, missing=mis.deldat))

     res
   } ## end function add.missing.deldat



  ## clean data
  cleandata.ReadData <- function(dat, datafile, marker.names, line.names)
  {
    ## Internal function
    ## Purpose: to identify delation, give a summary of the data read, and abbreviate 
    ##          the marker and line names
    ## Args:    dat  data frame 
    ##          datafile containing file name
    ## Outputs:  data matrix of type "delmap.data" with abbreviated marker & line names
    ##           and recoded genotypes to be 1,0, NA


  ## identify deletions and non-deletions
  ## (I have assumed deletions are fairly rare)
  genos <- as.vector(table(dat))
  nodelgeno <- names(table(dat))[which(genos== max(genos))]
  delgeno <- names(table(dat))[which(genos == min(genos))]

  cat("\n\n")
  cat(" Deletion data has been read in from the file:", datafile, "\n")
  cat(" The original data contains: \n")
  cat("        ",  nrow(dat), " rows (lines) \n")
  cat("        ",  ncol(dat), " columns (marker loci) \n")
  cat("        ",  round(sum(is.na(dat))/(nrow(dat)*ncol(dat)), 2), "% of the genotypes are missing\n")
  cat("         Deletions are assumed to be of genotype ", delgeno, "\n")
  cat("         Non-deletions are assumed to be of genotype", nodelgeno, "\n")

  ## recode deletion data to be 0/1 data
  recoded.dat <- matrix(data=0, nrow=nrow(dat), ncol=ncol(dat))
  recoded.dat[which(dat==delgeno, arr.ind=T)] <- 1
  recoded.dat[is.na(dat)] <- NA

  ## abbreviate marker names
  new.mrk.nms <- paste("M", 1:ncol(dat), sep="")
  colnames(recoded.dat) <- new.mrk.nms
  if(marker.names) colnames(recoded.dat) <- colnames(dat)


  ## abbreviate line names
  new.line.nms <- paste("L", 1:nrow(dat), sep="")
  rownames(recoded.dat) <- new.line.nms
  if(line.names) rownames(recoded.dat) <- rownames(dat)

  # new class of object where object contains matrix data with special attributes
  class(recoded.dat) <- "delmap.data"

  return(recoded.dat)
  }  ## end function cleandata.ReadData






####----------------------------------##
##   EXTERNAL FUNCTION DECLARATION  ##
###-----------------------------------##



SimDeletions <- function(numlines=NULL, plines=0.5, nummarkers=NULL, labels=NULL, Enumdel=4, 
                         p.missing=NULL, seed=NULL)
{

  ## Purpose:        to simulate deletion data.
  ## Functionality:  Generates a 2-dimensional array. The rows are the plant lines and 
  ##                 the columns are the genetic marker loci. The array contains marker 
  ##                 gentoypes 1/0 where 1 denotes a deletion (absent) and 0 denotes 
  ##                 a non-deletion (present).  
  ##
  ##
  ## Args
  ## numlines   -  number of HIB plants in population
  ## plines     -  prob of a line carrying a deletion
  ## nummarkers -  number of marker loci 
  ## labels     -  marker names (optional)
  ## Enumdel    -  expected number of deletions for poisson distribution of deletions
  ## p.missing  -  prob of a genotype being missing
  ## seed       -  a numeric value to initialize the pseudo-random number generator
  ##
  ## Returns    -  returns two delmap.data objects, one with all the genotypes observed
  ##               and the other with missing genotypes (if p.missing is not null). 
 

  ## set seed if not null
  if(!is.null(seed)) set.seed(seed)

  ## create 2D data structure
  deldat <- create.structure.SimDeletion(numlines, nummarkers, labels) 



   ## generate deletions by
   ##  1. sampling lines that contain deletions
   ##  2. sampling length of deletions
   ##  3. sampling location of deletions. 

  
  res <- step1.simulation(numlines, plines)

  deldat <- with(res, step23.simulation(rnumlines, Enumdel, numlines, nummarkers, chosen.lines, deldat ))

  
   return(add.missing.deldat(p.missing, deldat))

}  ## end function  SimDeletions

print.delmap.data <- function(x, ...)
{
  cat(" A summary of this object is: \n\n")
  cat(" The (recoded) deletion data are: \n")
  write.it(x)
  write.mrknames(x)
  write.linenames(x)
} ## end function print.delmap.data

ReadData <- function(datafile, line.names=FALSE,na.strings="NA", marker.names=FALSE, csv=FALSE)
{
 ## Purpose:        to read in deletion data.
 ## Functionality:  Reads in space and comma seperated data where marker/line names
 ##                 may or maynot be present. 

 if(missing(datafile))
   stop(" Input file cannot be missing.")

  ## check inputs
  checkinputs.ReadData(datafile, line.names, marker.names)

 fsep <- "";  if(csv)  fsep <- ","
 col.number <- NULL; if(line.names) col.number <- 1

  ## read data from file
  dat <- read.data(file=datafile, header= marker.names, sep= fsep, na.strings=na.strings, row.names=col.number)

  ## clean data (abbreviated marker & line names, of type delmap.data, and genotypes are
  ##    0, 1, NA
  cleandata.ReadData(dat, datafile, marker.names, line.names)
} ## end ReadData






