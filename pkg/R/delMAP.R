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
##   2. sort out missing genotypes and CleanData - ignore ?


  
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



UniqueCols <- function(mat) {
   ## Internal function
   ## Purpose: to identify unique columns
   ##          Returns only unique columns and a list of those marker loci removed. 
   ## mat:     a data matrix or delmap.data object



    # ucols will store the unique columns
    cols <- list()
    ncols <- dim(mat)[2]
    ucols <- 1:ncols

    # Iterate through columns comparing it with all columns
    # except ones previously checked.
    for(i in 1:(ncols-1)) {

        # Use to keep track of collapsed markers
        cols[[colnames(mat)[i]]] <- 0
        k <- 1

        # Iterate through each marker, if not unique, assign NA
        # to be removed at return stage.
        for(j in (i+1):ncols) {
            if(!is.na(ucols[j])) { 
                if(all(mat[,i]==mat[,j])) {
                    ucols[j] <- NA
                    cols[[colnames(mat)[i]]][k] <- colnames(mat)[j]
                    k <- k+1
                }
            }
        }

        if(length(cols[[colnames(mat)[i]]]) == 1) cols[[colnames(mat)[i]]] <- NULL
    }
    return(list(mat = mat[,ucols[!is.na(ucols)]], cols = cols))
}





  write.mrknames <- function(x)
  {
    ## Internal function
    ## Purpose: to write out marker names
    cat("\n")
    cat(" Marker names (in column data order) are: \n")
    if(length(rownames(x))  < 9)
        cat(" ", rownames(x),      "\n")
      else
        cat(" ", rownames(x)[1:9], "... \n")

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



CreateDistMatrix <- function(dmat=NULL)
{
## Purpose:  Use Manhattan distance measure on marker loci (columns)
## Args:     matrix or delmap.data object
## Return:   distance object for deletion data

   if(is.delmap(dmat)){
     class(dmat) <- "matrix"
   } else {
    if(is.matrix(dmat))
    {
      dmat <- as.delmap(dmat)
    } else stop(" Object must be of type matrix or delmap.data.")
   }
   if(any(is.na(dmat))) stop("NA values not allowed when calculating distances.")


    return(dist(t(dmat), method="manhattan"))

} ## end function CreateDistMatrix


CleanData <- function(dmat=NULL, ignoreNA=TRUE)
{
  ## Purpose:   CleanData 
  ##            Rows and columns without deletions are removed
  #             Default is to remove rows with only a single NA. 
  #             but this option is turned off by setting ignoreNA=FALSE

    if (!is.delmap(dmat)) 
        stop("Error! object not of class delmap.data")
    indx <- which(is.na(dmat), arr.ind = TRUE)


    if(ignoreNA)
    { ## ignoring NA's 
     indx <- which(rowSums(dmat, na.rm=TRUE) == 0) 
     if(length(indx) > 0) dmat <- dmat[-indx,]
     indx <- which(colSums(dmat, na.rm=TRUE)==0)
     if(length(indx) > 0) dmat <- dmat[, -indx]
    }  else {
      ## taking NA's into account in the count
     indx <- which(rowSums(dmat, na.rm=TRUE) == 0 & rowSums(is.na(dmat)) <2)
     if(length(indx) > 0) dmat <- dmat[-indx,]
     indx <- which(colSums(dmat, na.rm=TRUE) == 0 & colSums(is.na(dmat)) <2)
     if(length(indx) > 0) dmat <- dmat[,-indx]
    }     

    class(dmat) <- "delmap.data"
    dmat 

}  ## CleanData



PermuteCols <- function(dmat=NULL)
{
 if(class(dmat) != "delmap.data") stop("Error! object not of class delmap.data")

  cl <- class(dmat)
  ## permute columns
  indx <- sample(1:ncol(dmat), ncol(dmat), replace=F)
  a <- dmat[, indx]
  class(a) <- cl
  a
}



IdentifyMarkerBlocks <- function(dmat = NULL)
{
  ## Purpose:  to identify separate blocks of markers. 
  ##           In doing this, I am assuming that the data is 
  #            in the true ordering. Some patterns of deletion data
  ##           can lead to independent blocks of markers being formed.
  ##           These blocks can be in different orders. 
  ## Args:     dmat -  matrix of deletion data in true marker order
  ## Returns:  integer vector of blocks


  marker.groupings <- rep(FALSE, ncol(dmat)-1)
  ## Determining separate groupings of markers
  for(ii in 2:ncol(dmat))
    if (dmat[,ii-1] %*% dmat[, ii] > 0)
        marker.groupings[ii-1] <-  TRUE

  indx <- 1
  group.indx <- rep(NA, ncol(dmat))
  group.indx[1] <- indx
  for(ii in 2:ncol(dmat))
  {
    if(!marker.groupings[ii-1])
        indx <- indx + 1
    group.indx[ii] <- indx
  }
  return(group.indx)
} ## end function

IdentifyMarkerOrd <- function(dmat = NULL)
{

  ## Purpose:  to identify all candidate marker orderings
  ##           for each marker block of loci. 
  ## Args:     dmat  matrix of marker deletion data (0,1)
  ##           ord  ordering from seriate function
  ## Returns:  list where each element is a integer vector or possible orderings


  ## Order rows/lines based on deletion pattern
  d <- dmat
  class(d) <- "matrix" 
  sord <- seriate(d, method="BEA")

  ## identify any marker blocks
  blocks <- IdentifyMarkerBlocks(dmat)

  list.orders <- list(blocks=blocks, orders=NULL)

  for( jj in unique(blocks) )
  {

     indx <- which(blocks==jj)
     if(length(indx)==1)
     {
      list.orders$orders[[jj]] <- 1
     } else {
     dat <- dmat[, indx]
     ## determine: change in number of deletions across loci
     num.deletions <- colSums( dat[get_order(sord[1]),]   )
     change.in.number.of.deletions  <- rep(FALSE, ncol(dat)-1)
     for(ii in 2:ncol(dat))
        if( abs(num.deletions[ii] - num.deletions[ii-1]) > 0)
            change.in.number.of.deletions[ii-1] <- TRUE

     ## determine: change in lines with deletions
     change.in.lines.with.deletions <- rep(FALSE, ncol(dat)-1)
     for(ii in 2:ncol(dat))
          if( any((dat[,ii-1] + dat[,ii])==1) )
              change.in.lines.with.deletions[ii-1] <- TRUE

     ## if 0, then false false and marker not uniquely positioned
     col.sums <-  colSums(  matrix( data= c(change.in.number.of.deletions,
                          change.in.lines.with.deletions), byrow=T, nrow=2) )

     ## identify markers that cannot be uniquely positioned. 
     indx <- 1
     orderings <- rep(0, ncol(dat))
     orderings[1] <- indx
     for(ii in 2:ncol(dat))
     {
       if( col.sums[ii-1] != 0)
            indx <- indx + 1
       orderings[ii] <- indx
     }
     list.orders$orders[[jj]] <- orderings
     } ## end if else
    } ## end jj
    return(list.orders)

} ## end function IdentifyMarkerOrd




ImputeMissingGeno <- function(dmat=NULL)
{
  ## Assign 0,1 values to missing genotypes

 ## identify elements that have missing genotypes (NA)
 indx <- which(is.na(dmat), arr.ind=TRUE)
 dmat[indx] <- sample(unique(dmat)[!is.na(unique(dmat))], nrow(indx), replace=T)
 dmat
}  ## end function ImputeMissingGeno



plot.delmap.data <- function(x, main=NULL,...)
{


    class(x) <- "matrix"
    dimc <- ncol(x)
    dimr <- nrow(x)
    plot(x, xlim = c(1, dimc), ylim = rev(c(1, dimr)), xlab = "", 
        ylab = "", type = "n", axes = FALSE, frame = TRUE, bty="n", oma=c(1,1,1,1))
    if (!is.null(main)) 
        title(main = main, line = 3)
    axis(3, 1:dimc, dimnames(x)[[2]], cex.axis = 0.5, las = 3, mgp=c(0,1,0))
    axis(2, 1:dimr, dimnames(x)[[1]], cex.axis = 0.5, las = 1, mgp=c(0,0.5,-0.5))
    segments(x0 = 0.5:(dimc + 0.5), y0 = 0.5, y1 = max(dimr) + 0.5, col="grey")
    segments(x0 = 0.5, y0 = 0.5:(dimr + 0.5), x1 = max(dimc) + 0.5, col="grey")
    indx <- which(is.na(x), arr.ind = TRUE)
    points(indx[, 2], indx[, 1], cex=1)
    indx <- which(x == 1, arr.ind = TRUE)
    points(indx[, 2], indx[, 1], pch = 19, cex=1)
    t <- c(0, findInterval(unique(IdentifyMarkerBlocks(x)), IdentifyMarkerBlocks(x)))
    segments(x0 = t + 0.5, y0 = 0.5, y1 = max(dimr) + 0.5)
    axis(1, (t[-length(t)] + t[-1])/2, sprintf("B%d", unique(IdentifyMarkerBlocks(x))), 
        las = 3, tick = FALSE, mgp=c(0,0,0),cex.axis = 0.75)



#    class(x) <- "matrix"    
#    dimc <- ncol(x)
#    dimr <- nrow(x)
#    plot(x, xlim=c(1,dimc), ylim=rev(c(1,dimr)), xlab="", ylab="", type="n", axes=FALSE, frame=TRUE)
#    if(!is.null(main)) title(main=main, line=3)
#    axis(3, 1:dimc, dimnames(x)[[2]], cex.axis = 0.7, las=3)
#    axis(2, 1:dimr, dimnames(x)[[1]], cex.axis = 0.7, las=1)
#    # set grid lines
#    abline(v=0.5:(dimc+0.5), col="grey")
#    abline(h=0.5:(dimr+0.5), col="grey")
#    # add data
#    indx <- which(is.na(x), arr.ind=TRUE)
#    points(indx[,2], indx[,1])
#    indx <- which(x==1, arr.ind=TRUE)
#    points(indx[,2], indx[,1], pch=19)
#    # determine and plot blocks
#    t <- c(0, findInterval(unique(IdentifyMarkerBlocks(x)),IdentifyMarkerBlocks(x)))
#    abline(v=t+0.5)
#    axis(1, (t[-length(t)]+t[-1])/2, sprintf("B%d", unique(IdentifyMarkerBlocks(x))), las=3, tick=FALSE)

#     }
#   }

} ## end function plot.delmap.data

as.delmap.matrix <- function(x)
{
  # convert matrix object into delmap.data object.
  # matrix is given generic marker and line names. 
  
  # checks
  if(ncol(x) < 2) stop("Matrix must contain more than a single column.")
  if(nrow(x) < 2) stop("Matrix must contain more than a single row.")
  if(length(unique(as.vector(x))[!is.na(unique(as.vector(x)))]) > 2) stop("Matrix contains more than two unique values.")

  # convert matrix scores into 0/1 where we assume deletions are least prevalent. 
  mode(x) <- "factor"  
  mode(x) <- "numeric" ## trick to give me numbers
  t <- table(x)
  if(t[1] > t[2] )
  {
     x[which(x== as.numeric(names(t[1])), arr.ind=TRUE)] <- 0
     x[which(x== as.numeric(names(t[2])), arr.ind=TRUE)] <- 1
  } else {
     x[which(x== as.numeric(names(t[2])), arr.ind=TRUE)] <- 0
     x[which(x== as.numeric(names(t[1])), arr.ind=TRUE)] <- 1
 } 
 
  ## create generic marker names
  colnames(x) <- paste("M", 1:ncol(x), sep="")

  ## abbreviate line names
  rownames(x)  <- paste("L", 1:nrow(x), sep="")

  # new class of object where object contains matrix data with special attributes
  class(x) <- "delmap.data"
  
  x
}

is.delmap <- function(dmat) inherits(dmat, "delmap.data")

as.delmap <- function(x) 
  UseMethod("as.delmap")


as.delmap.default <- function(x)
{
  if(is.delmap(x))
     x
  else
    {
      stop(" Only objects of type matrix can be converted into delmap.data objects.")
      NULL
    } 

}

#===================================================
# Function to collapse like markers into one
# marker column. If a column contains missing values
# i.e. NA then it is not collapsed because don't know
# if NA is 0 or 1.
# Returns a list with the collapsed data and info on which
# columns were collapsed
# dmat: a delmap data object
#===================================================
Collapse <- function(dmat=NULL) {

   class(dmat) <- "matrix"

   # The current ordering of dmat
   orig.ord <- 1:ncol(dmat)
   names(orig.ord) <- colnames(dmat)

   indx <- which(is.na(dmat), arr.ind=TRUE)

   # Assume no missing data
   miss <- NULL
   nomiss <- dmat

   # If a marker has missing data, assign to miss,
   # otherwise assign to no miss.
   if(length(indx)!=0) {
       miss <- dmat[,unique(indx[,2])] 
       nomiss <- dmat[,-unique(indx[,2])]
   }

   uniq <- UniqueCols(nomiss)

   # Combine the unique columns with the
   # miss matrix
   ret <- cbind(uniq$mat, miss)
   ret <- ret[,names(sort(orig.ord[colnames(ret)]))]
   class(ret) <- "delmap.data"
   return(list(data = ret, collapsed = uniq$cols))
}

