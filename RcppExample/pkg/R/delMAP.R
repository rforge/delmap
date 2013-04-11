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

#===================================================
# This function moves the configuration to a neighbour
# in the configuration space
# x: the data vector
# u: the elements of x that are to be swapped
#===================================================
Perturb <- function(x, u) {

    # Turn 0's to 1's and 1's to 0's
    for(i in 1:length(u)) {
        if(x[u[i]] == 0 ) 
            x[u[i]] <- 1
        else
            x[u[i]] <- 0
    }
    return(x)
}



#===================================================
# Sort randomised columns when there are no missing
# data values, could be combined into the Order function
# but have kept it seperate for now
# dmat :  delmap data
# option: 1 for Manhattan dist, otherwise delmap dist
# n: number of iterations (due to it being heuristic)
#===================================================
Sort <- function(dmat = NULL, n = 1000)
{


    cbest <- CleanData(dmat)
    dbest <- CreateDistMatrix(cbest[,])
    lbest <- criterion(dbest, method="Path_length")[1]
    obest <- seriate(dbest, method="TSP", control=list(method="2-opt"))
    ordbest <- cbest[,get_order(obest)]

    # used to keep track of how lengths are changing
    lengths <- rep(NA, n)

    # Ordering obtained by Seriate TSP is not necessarily the same each time
    # Iterate through to make sure we get the best, 1000 probably more than
    # needed
    for(i in 1:n) {
        d <- CreateDistMatrix(ordbest[,])
        o <- seriate(d, method="TSP", control=list(method="2-opt"))
        ord <- ordbest[,get_order(o)]
        l <- criterion(CreateDistMatrix(ord[,]), method="Path_length")[1]

        if(l < lbest)
        {
            lbest <- l
            ordbest <- ord
        }
        lengths[i] <- l
    }
    return(list(ord=as.delmap(ordbest), len=lengths))
}



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
    if(length(colnames(x))  < 9)
        cat(" ", colnames(x),      "\n")
      else
        cat(" ", colnames(x)[1:9], "... \n")

  } ## end function write.mrknames


  write.linenames <- function(x)
  {
    ## Internal function
    ## Purpose:  write out line names 

    cat("\n")
    cat(" Line names (in row data order) are: \n")
    if(length(rownames(x)) < 5)
       cat(" ", rownames(x),  "\n\n")
     else
       cat(" ", rownames(x)[1:5], "... \n\n")

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


CleanData <- function(dmat=NULL, ignoreNA=FALSE)
{
  ## Purpose:   CleanData 
  ##            Removing of noninformative rows/columns.
  ##            If ignoreNA=TRUE, rows and columns are removed without regard to NA's if they carry no deletions
  ##            IF ignoreNA=FALSE, only rows without deletions and/or columns with no NA's and 
  ##            no deletions are removed. That is, a column is only removed if it carrys no deletions and has 
  ##            no NA's 
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
     indx <- which(rowSums(dmat, na.rm=TRUE) == 0 )
     if(length(indx) > 0) dmat <- dmat[-indx,]
     indx <- which(colSums(dmat, na.rm=TRUE) == 0 & colSums(is.na(dmat)) ==0)
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
  ##
  ## Note:     All columns with no deletions are assigned to be frome the same block
  ##           NA's not allowed


  if(any(is.na(dmat))) stop("Error.  NA's not allowed.")
  if(!is.delmap(dmat)) stop("Error.  Object must be of class delmap.data.")

  marker.groupings <- rep(FALSE, ncol(dmat)-1)
  group.indx <- rep(NA, ncol(dmat))
  zeroblockindx <- NULL
 
  ## Determining separate groupings of markers
  for(ii in 2:ncol(dmat))
    if (dmat[,ii-1] %*% dmat[, ii] > 0)
        marker.groupings[ii-1] <-  TRUE

  if(sum(dmat[,1])==0) { 
     zeroblockindx <- 1
     blockindx <- 1
     group.indx[1] <- zeroblockindx
   } else {
     blockindx <- 1
     group.indx[1] <- blockindx
   }
  for(ii in 2:ncol(dmat))
  {
    if(sum(dmat[,ii])==0){
       if(is.null(zeroblockindx)) {
         # zeroblock indx not set yet
         zeroblockindx <- blockindx + 1
         blockindx <- blockindx + 1
        }
       group.indx[ii] <- zeroblockindx
     } else {
      if(!marker.groupings[ii-1])
          blockindx <- blockindx + 1
      group.indx[ii] <- blockindx
    }
  }
  names(group.indx) <- colnames(dmat)
  return(group.indx)
} ## end function

IdentifyMarkerOrd <- function(dmat = NULL)
{

  ## Purpose:  to identify all candidate marker orderings
  ##           for each marker block of loci. 
  ## Args:     dmat  matrix of marker deletion data (0,1)
  ##           ord  ordering from seriate function
  ## Returns:  list where each element is a integer vector or possible orderings


  if(!is.delmap(dmat)) stop("Error. Object must be of class delmap.data.")


  ## Order rows/lines based on deletion pattern
  d <- dmat
  class(d) <- "matrix" 
  sord <- seriate(d, method="BEA")

  ## identify any marker blocks
  blocks <- IdentifyMarkerBlocks(dmat)

  list.orders <- list(blocks=blocks, orders=list())

  for( jj in unique(blocks) )
  {

     cindx <- which(blocks==jj)
     if(length(cindx)==1)
     {
      list.orders$orders[[jj]] <- 1
      names(list.orders$orders[[jj]]) <- names(blocks)[cindx]
     } else {
     dat <- dmat[, cindx]
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
     names(list.orders$orders[[jj]]) <- colnames(dmat)[cindx]
     } ## end if else
    } ## end jj
    return(list.orders)

} ## end function IdentifyMarkerOrd




ImputeMissingGeno <- function(dmat=NULL,  uniform=TRUE)
{
   if(!is.delmap(dmat)) stop("Object must be of class delmap.data.")

   indx <- which(is.na(dmat), arr.ind=TRUE)
   if(uniform)
      dmat[indx] <- sample(unique(dmat)[!is.na(unique(dmat))], nrow(indx), replace=T)

   if(!uniform) {
     
      for(ii in 1:nrow(indx))
      {
        ro <- indx[ii,1]; co <- indx[ii,2]
        if (co > 1 & co < ncol(dmat) ) {
           jjr <- jjl <- co
           
           while(is.na(dmat[ro,jjr])  & jjr != ncol(dmat)) jjr <- jjr+1
           while(is.na(dmat[ro,jjl]) &  jjl != 1)           jjl <- jjl- 1
           if (sum(dmat[ro, c(jjl, jjr)], na.rm=TRUE) == 0)  dmat[ro,co] <- sample(0:1, 1, prob=c(0.95, 0.05))
           if (sum(dmat[ro, c(jjl, jjr)], na.rm=TRUE) == 1)  dmat[ro,co] <- sample(0:1, 1, prob=c(0.5, 0.5))
           if (sum(dmat[ro, c(jjl, jjr)], na.rm=TRUE) == 2)  dmat[ro,co] <- sample(0:1, 1, prob=c(0.05, 0.95))
         } ## end if

        if (co==1) {
           jjr <- co
           while(is.na(dmat[ro,jjr])  & jjr != ncol(dmat)) jjr <- jjr+1
           if (dmat[ro, jjr] == 0)  dmat[ro,co] <- sample(0:1, 1, prob=c(0.95, 0.05))
           if (dmat[ro, jjr] == 1)  dmat[ro,co] <- sample(0:1, 1, prob=c(0.5, 0.5))
         } # end if
         if (co == ncol(dmat)) {
           jjl <- co
           while(is.na(dmat[ro,jjl]) & jjl != 1) jjl <- jjl - 1
           if (dmat[ro, jjl] == 0)  dmat[ro,co] <- sample(0:1, 1, prob=c(0.95, 0.05))
           if (dmat[ro, jjl] == 1)  dmat[ro,co] <- sample(0:1, 1, prob=c(0.5, 0.5))
         } # end if
      } # end for
   }  ## end if !uniform
   return(dmat)

}  ## end function 








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
    t <- c(0, findInterval(unique(IdentifyMarkerBlocks(as.delmap(x))), IdentifyMarkerBlocks(as.delmap(x))))
    segments(x0 = t + 0.5, y0 = 0.5, y1 = max(dimr) + 0.5)
    axis(1, (t[-length(t)] + t[-1])/2, sprintf("B%d", unique(IdentifyMarkerBlocks(as.delmap(x)))), 
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
  # matrix is given generic marker and line names if col and row names are null.
  
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
  if(is.null(colnames(x))) colnames(x) <- paste("M", 1:ncol(x), sep="")

  ## abbreviate line names
  if(is.null(rownames(x))) rownames(x)  <- paste("L", 1:nrow(x), sep="")

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



#===================================================
# Function to find an approximately optimal ordering
# according to Hamiltonion path length.
#
# dmat:   the delmap data
# n:      the no. of iterations to perform
# option: 1 for Manhattan distance, otherwise, delmap
# w:      the "window" to work with i.e how many values
#         change during each iteration
# T:      value of temperature
#===================================================
Order <- function(dmat=NULL, n=1000,  w=1, T=1000) {
    if (class(dmat) != "delmap.data")
        stop("Error! object not of class delmap.data")

    # Record original ordering in case we want it
    orig.ord <- 1:ncol(dmat)
    names(orig.ord) <- colnames(dmat)

    # identify missing values
    indx <- which(is.na(dmat[,]), arr.ind=TRUE)

    # If no missing values just use sort
    if(length(indx) < 1) {
        print("No missing values! Used Sort() instead.") 
        res <- Sort(dmat)
        return(res)
    }
    
    # used to store locations of missing values          
    ind.na <- list(row=rownames(dmat)[indx[,1]], col=colnames(dmat)[indx[,2]], value=dmat[indx])

    # Set initial configuration and store
    S <- ImputeMissingGeno(dmat)
    Sc <- CleanData(S)
    Sd <- CreateDistMatrix(Sc)
    Sl <- criterion(Sd, method="Path_length")[1]
    So <- seriate(Sd, method="TSP", control=list(method="2-opt"))
    Sord <- Sc[,get_order(So)]

    # Stores the imputed values of missing data
    bestvals <- sapply(FUN = function(x) S[ind.na[[1]][x], ind.na[[2]][x]], 1:length(ind.na[[1]]))
    ind.na$value <- bestvals

    # Set arbitrarily as best current ordering
    Sbest <- S
    Sbestc <- Sc
    Sbestd <- Sd
    Sbestl <- Sl
    Sbesto <- So
    Sbestord <- Sord

    # Used to store new length, accepted lengths and best length at each iteration
    lengths <- rep(NA,n)
    lengths.a <- rep(NA,n)
    Sbestl.i <- rep(NA, n)

    # Iterate through algorithm n times
    for(i in 1:n) {

        # determine neighbour to investigate
        currvals <- sapply(FUN = function(x) S[ind.na[[1]][x], ind.na[[2]][x]], 1:length(ind.na[[1]]))
        u <- sample(length(currvals), size=w, replace=FALSE)

        # Move to neighbour in configuration space
        newvals <- Perturb(currvals, u)
        Snew <- S
        for(k in 1:length(u)) {
            Snew[ind.na[[1]][u[k]], ind.na[[2]][u[k]]] <- newvals[u[k]]
        }

        # calculate "energy" for neighbour
        Snewc <- CleanData(Snew)
        Snewd <- CreateDistMatrix(Snewc)
        Snewo <- seriate(Snewd, method="TSP", control=list(method="2-opt"))
        Sneword <- Snewc[,get_order(Snewo)]
        Snewdord <- CreateDistMatrix(Sneword)
        Snewl <- criterion(Snewdord, method="Path_length")[1]
        lengths[i] <- Snewl

        # Make the comparisons and accept new
        # configuration with probability P(Sl, Snewl, T).
        v <- runif(1,0,1)
        if(exp((Sl - Snewl)/T) > v) {
            S <- Snew
            Sc <- Snewc
            Sd <- Snewd
            Sl <- Snewl
            So <- Snewo
            Sord <- Sneword
            lengths.a[i] <- Snewl
        }
        else
            lengths.a[i] <- 0

        # If new configuration is better than current
        # best, save as new best configuration.
        if(Snewl < Sbestl) {
            Sbest <- Snew
            Sbestc <- Snewc
            Sbestd <- Snewd
            Sbestl <- Snewl
            Sbesto <- Snewo
            Sbestord <- Sneword
            ind.na$value <- sapply(FUN = function(x) Sbest[ind.na[[1]][x], ind.na[[2]][x]], 1:length(ind.na[[1]]))

            # Plot to get an idea of how things are going
            plot(as.delmap(Sbestord), main=sprintf("Plot at iteration %d", i))
        }

        Sbestl.i[i] <- Sbestl

        # Update parameters
        T <- T*0.99
        i <- i + 1
    }

    class(Sbestc) <- "delmap.data"
    class(Sbestord) <- "delmap.data"

    plot(Sbestord, main="Final Plot")

    # Return the unordered and ordered data along with lengths and best lengths
    list(unord=Sbestc[,intersect(names(orig.ord), names(get_order(Sbesto)))], 
         res=Sbestord, ord=Sbesto, index=ind.na, lengths=lengths, acc=lengths.a, best=Sbestl.i)
}


#DeletionMappingOLD <- function(dmap=NULL, niterates=100, nwithin=100,temp=1000, psampled=0.1, method="concorde",...)
#{
#
#    if(!is.delmap(dmap)) stop("Object must be of class delmap.data.")
#   
#    ##-----------------------------##
#    ##  Initialization             ##
#    ##-----------------------------##
#    dmap <- CleanData(dmap, ignoreNA=FALSE)
#    odmap <- ImputeMissingGeno(dmap)
#    o.order <- 1:ncol(dmap)
#    names(o.order) <- colnames(dmap)
#    i.missing <- which(is.na(dmap), arr.ind=TRUE)  # index of missing values
#    v.missing <- odmap[i.missing]                  # values assigned to missing
#    otlength <- 100000   # tour length set to arbritarily larger value to ensure not accepted
#
#    bestmap <- odmap
#    besttlength <- otlength
#
#
#    for (ii in 1:niterates)
#    {
#       
#      for(jj in 1:nwithin)
#      {
#       # create new realization
#         # identify missing values to replace with new sampled values
#         rindx <- sample(1:nrow(i.missing), round(nrow(i.missing)*psampled), replace=FALSE)
#         nvals <- sample(unique(dmap[!is.na(dmap)]),length(rindx), replace=TRUE)
#         # form new realization from old dmap 
#         ndmap <- odmap  # new dmap
#         ndmap[matrix(i.missing[rindx,],ncol=2,byrow=FALSE) ] <- nvals 
#         class(ndmap) <- "delmap.data"
#
#       # find best ordering of new realization and its touring length
#        D <- CreateDistMatrix(ndmap)
#        tD <- as.TSP(D)  # in TSP format
#        tD <- insert_dummy(tD, label="cut")  # adding extra dummy city
#        solu   <- solve_TSP(tD, method=method)
#        ntlength <- attr(solu, "tour_length")
#        n.order <- cut_tour(solu, "cut")  # new marker ordering
#
#       # test if we should move to new realization
#       rnd <- runif(1,0,1)
#       MHprob <- exp( -1*(ntlength-otlength)/temp)
#       if (MHprob >= 1 | rnd < MHprob)
#       { # accept new realization
#         # trick part - need to rearrange i.missing and v.missing to reflect to column ordering
#         lookuptab <- data.frame(old=1:ncol(odmap), new=n.order)
#         indx <- with(lookuptab, match( i.missing[,2], old))         
#         i.missing[,2] <- lookuptab$new[indx]
#
#         ndmap <- odmap[,n.order]
#         class(ndmap) <- "delmap.data"
#         odmap <- ndmap
#         o.order <- n.order
#         otlength <- ntlength
#    
#         #keep track of best realization
#         if(ntlength < besttlength)
#         {
#           bestmap <- ndmap
#           besttlength <- ntlength
#         }
#
#       } # end if 
#      } # end for jj nwithin 
#       # geometric update temperature
#       temp <- 0.95 * temp
#
#     #  plot(bestmap)
#      } # end for ii niterates
#         list(dmap=bestmap, tlength=besttlength)
#}  # end function

iDeletionMapping <- function(dmap,  psampled, odmap, temp, method, otlength, besttlength, bestmap, blockstr, mrkord, bestdist)
{
       # impute missing marker genotypes
       ndmap <- ImputeMissingGeno(dmat=dmap, uniform=FALSE)


       # find best ordering of new realization and its touring length
       D <- CreateDistMatrix(ndmap)
       tmpD <- as.matrix(D)
       indxr <- which(rownames(tmpD)=="cut")
       indxc <- which(colnames(tmpD)=="cut")
       tmpD[indxr,] <- 0
       tmpD[,indxc] <- 0
       D <- as.dist(tmpD)

       tD <- as.TSP(D)  # in TSP format
       ### tD <- insert_dummy(tD, label="cut")  # adding extra dummy city
       solu   <- solve_TSP(tD, method=method)
       ntlength <- attr(solu, "tour_length")
       n.order <- cut_tour(solu, "cut")  # new marker ordering
       n.order <- c(n.order, ncol(ndmap)) # adding "cut" back in
       # test if we should move to new realization
       rnd <- runif(1,0,1)
       MHprob <- exp( -1*(ntlength-otlength)/temp)
       if (MHprob >= 1 | rnd < MHprob)
      { # accept new realization
        # rearrange dmap to be new order
         dmap <- dmap[, n.order]  # reorder data with missing data
         ndmap <- ndmap[,n.order] # reorder  newly imputed data
         class(dmap) <- class(ndmap) <- "delmap.data"

         odmap <- ndmap  # assign newly imputed data to current data
         otlength <- ntlength
        #keep track of best realization
         if(ntlength < besttlength)
         {
           bestmap <- ndmap
           besttlength <- ntlength
           blockstrbest <- IdentifyMarkerBlocks(bestmap)
           bestdist <- RearrangeDist(blockstr, mrkord, blockstrbest)
         }
        } ## end if

        if(!is.null(blockstr) & !is.null(mrkord))

        blockstrobs <- IdentifyMarkerBlocks(odmap)
        tmp <- RearrangeDist(blockstr, mrkord, blockstrobs) ## calculate distance metric
        print(tmp)

        res <- list(bestmap=bestmap, besttlength=besttlength, odmap=odmap, 
                     otlength=otlength, dmap=dmap, bestdist=bestdist)



        return(res)
}



RcppDeletionMapping <- function(dmap,  psampled, odmap, cooling, method, otlength, besttlength, bestmap, 
                                 niterates, nwithin, blockstr, mrkord, bestdist){
    .Call( "RcppDeletionMapping", dmap,  psampled,  odmap, cooling, method, 
                                  otlength, besttlength, bestmap, 
                                  niterates, nwithin, blockstr, mrkord, bestdist, PACKAGE="delmap")
}

DeletionMapping <- function(dmap=NULL, niterates=100, nwithin=100,cooling=0.99 , psampled=0.1, 
                            method="concorde", blockstr=NULL, mrkord=NULL, ...)
{
    if(!is.delmap(dmap)) stop("Object must be of class delmap.data.")
     if( cooling < 0 | cooling > 1) stop(" Cooling constant must be between 0 and 1")
    ##-----------------------------##
    ##  Initialization             ##
    ##-----------------------------##
    dmap <- CleanData(dmap, ignoreNA=FALSE)

    ## add cut to dmap
    n <- colnames(dmap)
    dmap <- cbind(dmap, rep(0, nrow(dmap)))
    colnames(dmap) <- c(n , "cut")
    class(dmap) <- "delmap.data"

    odmap <- ImputeMissingGeno(dmap)
    bestmap <- odmap
    otlength <- besttlength <- 1000000  # starting tour length.
    bestdist <- list("ntrueblocks"=0, "nobsblocks"=0, "mindist"=0, "minwrong"=0)

    res <- RcppDeletionMapping(dmap,  psampled, odmap,cooling, method, 
                             otlength, besttlength, bestmap, niterates, nwithin,
                             blockstr, mrkord, bestdist)

    return(list(bestmap=res$bestmap, besttlength=res$besttlength, t=res$t, bestdist=res$bestdist))

}  # end function


##--------------------------------##
## Under Construction             ##
##--------------------------------##

## workign on distance Metric

d <-   SimDeletions( numlines = 50, plines = 0.2, nummarkers = 50,
    Enumdel = 8, p.missing = 0.1,  seed = 1)

dc <- CleanData(d[["missing"]])


## obtain true blocking structure
mrks <- colnames(dc)
cindx <- match(mrks, colnames(d[["nomissing"]]))
dn <- as.delmap(d[["nomissing"]][,cindx])
trueblockstr <- IdentifyMarkerBlocks(dn)
truemrkord   <- IdentifyMarkerOrd(dn)$orders




#===================================================
# Calculates the L1 distance for a single symbol
# x: the symbol locations for the text string
# y: the symbol locations for the pattern string
#===================================================
l1Dist.single <- function(x, y) {
    sum(abs(x-y))
}

#===================================================
# Calculates the L1 distance of two strings
# Refer to Chapter 2 of pattern matching paper.
# x: the text string (resultant order from algorithm)
# y: the pattern string (known true ordering)
# i: attempt to account for strings of differing lengths
#===================================================
l1Dist <- function(x, y, i=1) {

    cost1 <- 0
    cost2 <- 0

    if(length(x) == length(y)) {

        # Extract index locations for each seperate symbol
        phix <- lapply(sort(unique(x)), function(k) which(x==k))
        phiy <- lapply(sort(unique(y)), function(k) which(y==k))
        phiz <- lapply(sort(unique(y)), function(k) which(rev(y)==k))

        # Calculate cost for both y and rev(y) i.e. z
        cost1 <- sum(sapply(1:length(phix), function(i) l1Dist.single(phix[[i]],phiy[[i]])))
        cost2 <- sum(sapply(1:length(phix), function(i) l1Dist.single(phix[[i]],phiz[[i]])))
    }

    # For what this is attempting to do: refer to 3.2 in "Pattern Matching
    # with address errors: rearrangement distances" The Case of Text and Pattern
    # of Different Sizes. Doesn't quite work how I want it to.
    else {

        # Extract index locations for each seperate symbol, but only
        # for the symbols that are in both strings
        phix <- lapply(sort(unique(x)), function(k) which(x==k))
        phiy <- lapply(sort(unique(x)), function(k) which(y[y%in%x]==k))
        phiz <- lapply(sort(unique(x)), function(k) which(rev(y[y%in%x])==k))
        
        
        y.i <- which(y==(y[y%in%x]))
        z.i <- which(rev(y)==(rev(y)[rev(y)%in%x]))

        # Calc distance for single symbol, but adjust the position index by i
        cost1 <- sum(sapply(1:length(phix), function(k) l1Dist.single(phix[[k]]+i-1,phiy[[k]]+y.i[1]-1)))
        cost2 <- sum(sapply(1:length(phix), function(k) l1Dist.single(phix[[k]]+i-i,phiz[[k]]+z.i[1]-1)))
    }
    return(min(cost1,cost2))
}



RearrangeDist <- function(blockstr=NULL, mrkord=NULL, blockstrobs=NULL)
{
## blockstr    vector object for true block structure (order specific)
## mrkord      list object of true marker ordering within block structure (order specific)
## blockstrobs vector object of realized block structure from simulated annealing  (order specific)

##-------------------------------------------------------------------------##
## Overview of Distance Metric                                             ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                           ##
##                                                                         ##
## Very challenging to invent a single univariate metric                   ##
## Have to cope with differences in block structure, differeces in         ##
## which markers appear in the blocks, the order of the markers within     ##
## a block, different numbers of markers within a block.                   ##
##                                                                         ##
## Example                                                                 ##
## Suppose the truth is as follows                                         ##
##  M1 M2 M3 M4 M5     M6 M7 M8 M9   -- marker names                       ##
##  1  1  2  2  3      1  1  2  2    -- pattern of identical markers       ##
##                                                                         ##
## Realized/observed ordering                                              ##
##  M1 M6 M5 M4   M2 M7    M3 M6 M8      M9                                ##
##   1  1  2  2   1  2     1  2  3       1                                 ##
##                                                                         ##
##                                                                         ##
## Multivariate distance metric                                            ##
## Dimensions  1. # of blocks                                              ##
##             2. total min rearrangement distance                         ##
##             3. number of wrongly blocked markers                        ##
##                                                                         ##
## # of blocks: (True)  2                                                  ##
##              (False) 4                                                  ##
##                                                                         ##
## total min rearrangement distance                                        ##
##  Block1:                                                                ##
##  compare M1  M4 M5                                                      ##
##          1   2  3                                                       ##
##  with                                                                   ##
##          M1  M5 M4                                                      ##
##          1   3  2                                                       ##
## Alot going on                                                           ##
## * first, for blocks 1, 2, 3, and 4, I identified the true block       ##
##      that it "best" matches.                                            ##
## * second, I am only considering those markers that are in common        ##
##      between the truth and the realized.                                ##
## * third, I am recoding the realized markers based on the recoding       ##
##      used in the truth (i.e. 1 3 2 for realized instead of              ##
##      1  2  2 for realized.                                              ##
##      this is because only the truth matters and I want to see how       ##
##      close my ordering is to the truth.                                 ##
##                                                                         ##
## minimum number of misplaced markers (i.e. in wrong blocks               ##
##                                                                         ##
##-------------------------------------------------------------------------##

##------------------------------------##
## Number of blocks                   ##
##------------------------------------##

## number of blocks - true and observed
res <- list()
res[["ntrueblocks"]] <- length(unique(blockstr))
res[["nobsblocks"]]  <- length(unique(blockstrobs))



##---------------------------------------------##
## Total minumum rearrangement distance        ##
## Total minimum wrongly placed markers        ##
##---------------------------------------------##


## calculate total rearrangement distance
total.dist <- 0  ## initalize total rearrangement distance
total.wrong <- 0 ## initalize total min number of wrongly placed markers
for(ii in unique(blockstrobs))
{
 ## identify which true block best matches the observed block (based on max number of matching markers)
 m.obs <- names(blockstrobs)[which(blockstrobs==ii)] 

 tb <- table(blockstr[match(m.obs, names(blockstr))])  ## frequency of markers in true blocks
 best.block <- as.numeric(names(which(max(tb)==tb)))  ## best block (may be more than one if ties

 vdist <- rep(NA, length(best.block))
 vwrong <- rep(NA, length(best.block))  ## wrongly placed markers
 for(jj in best.block)  ## allowing for ties of best block
 {
   ##  only keep matching (intersection of)  markers in m.obs and m
   m <- names(blockstr)[which(blockstr==jj)]
   m.obs.match <- m.obs[!is.na(match(m.obs, m))]
   m.match <- m[!is.na(match(m, m.obs))]
   ## recode common markers in m and m.obs but conditional on true marker recoding (mrkord)
   m.obs.recoded <- unlist(mrkord)[match(m.obs.match, names(unlist(mrkord)))]  ## using coding of true ordering
   m.recoded <- unlist(mrkord)[ match(m.match, names(unlist(mrkord)))]  ## using coding of true ordering
 
   vdist[which(best.block==jj)] <- l1Dist(m.obs.recoded, m.recoded)  ## rearrangement distance
   vwrong[which(best.block==jj)] <- length(m.obs) - length(m.obs.match)  ## number of wrongly placed markers
 }     ## end for jj
 total.dist <- total.dist + min(vdist)
 total.wrong <- total.wrong + min(vwrong)

 
} ## end for ii over unique observed blocks
res[["mindist"]] <- total.dist
res[["minwrong"]] <- total.wrong

return(res)

} ## end function







