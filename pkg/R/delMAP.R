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
##   2. sort out missing genotypes and cleandata - ignore ?


  
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
# dmat :  dmdat
# option: 1 for Manhattan dist, otherwise delmap dist
# n: number of iterations (due to it being heuristic)
#===================================================
Sort <- function(dmat = NULL, n = 1000)
{


    cbest <- cleandata(dmat)
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
    return(list(ord=as.dmdata(ordbest), len=lengths))
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
   ## mat:     a data matrix or dmdata object



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
   ## Purpose:     To create data structure for sim.dmdata
 
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
       if(nummarkers - delsize[ii] < 0) {
          ## entire row are 1's
          deldat[ chosen.lines[ii], ] <- 1
        } else {
         start.rindx <- sample( 1:(1+nummarkers- delsize[ii]), 1) 
         deldat[ chosen.lines[ii],  start.rindx:(start.rindx+delsize[ii] - 1)  ] <- 1
        } ## end if
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
      class(mis.deldat) <- "dmdatachrm"
      attr(mis.deldat, "imputed") <- which(is.na(mis.deldat), arr.ind=TRUE)
    }

    class(deldat) <- "dmdatachrm"
    attr(deldat, "imputed") <- NA
    ifelse( is.null(p.missing),
       res <- list(nomissing=deldat),
       res <- list(nomissing=deldat, missing=mis.deldat))

     res
   } ## end function add.missing.deldat



#  ## clean data
#  cleandata <- function(dat, file, names.pres )
#  {
#    ## Internal function
#    ## Purpose: to identify deletion, give a summary of the data read
#    ## Args:    dat  data frame 
#    ##          file containing file name
#    ##          names.pres logical for whether the labels are supplied or to be generated
#    ## Outputs:  data matrix of type "dmdata" with marker & line names
#    ##           and recoded genotypes to be 1,0, NA
#
#
#  genos <- as.vector(table(dat))
#  nodelgeno <- names(table(dat))[which(genos== max(genos))]
#  delgeno <- names(table(dat))[which(genos == min(genos))]
#
#  cat("\n\n")
#  cat(" Deletion data has been read in from the file:", file, "\n")
#  cat(" The original data contains: \n")
#  cat("        ",  nrow(dat), " rows (lines) \n")
#  cat("        ",  ncol(dat), " columns (marker loci) \n")
#  cat("        ",  round(sum(is.na(dat))/(nrow(dat)*ncol(dat)), 2), "% of the genotypes are missing\n")
#  cat("         Deletions are assumed to be of genotype ", delgeno, "\n")
#  cat("         Non-deletions are assumed to be of genotype", nodelgeno, "\n")
#
#  ## recode deletion data to be 0/1 data
#  recoded.dat <- matrix(data=0, nrow=nrow(dat), ncol=ncol(dat))
#  recoded.dat[which(dat==delgeno, arr.ind=T)] <- 1
#  recoded.dat[is.na(dat)] <- NA
#
#  return(recoded.dat)
#  }  ## end function 




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


  i.checks.subsetdata <- function(x=NULL, keep.chrm=NULL, drop.mrks=NULL, drop.rows=NULL, keep.mrks=NULL, keep.rows=NULL)
  {
    ## internal function to check arguements to subsetdata
    ## Args
    ##  x             an object of type dmdata
    ##  drop.mrks     a integer, numeric, or character vector of marker loci/columns to be dropped 
    ##  drop.rows     a integer, numeric, or character vector of rows to be dropped 
    ##  keep.mrks     a integer, numeric, or character vector of marker loci/columns to be kept 
    ##  keep.rows     a integer, numeric, or character vector of rows to be kept 

  ## checks
   ## check x
  if (!is.dmdata(x)  )
      stop("Object not of class dmdata")
  ## check chrm
  if(length(x) > 2) ## when there is more than a single chromosome
  {
     if(is.null(keep.chrm))
        stop(" At least one chromosome must be specified with keep.chrm.")
  }


  if(!is.null(keep.chrm)) {
      if(! (is.character(keep.chrm) | is.integer(keep.chrm) | is.numeric(keep.chrm) ) )
         stop(" chrm, when a vector, must be of class character, integer, or numeric. ")
  }

  ## check drop.mrks and drop.rows
  if(!is.null(drop.mrks) | !is.null(drop.rows)) {
     if(length(keep.chrm) > 1)
        stop(" When drop.mrks and/or drop.rows are specificed, only a single chromosome can be specified.")
  }

  ## check drop.mrks
  if(!is.null(drop.mrks)) {
       if(!(is.character(drop.mrks) | is.integer(drop.mrks) | is.numeric(drop.mrks)) )
         stop(" drop.mrks, when a vector, must be of class character, integer, or numeric. ")
  }

  ## check drop.rows
  if(!is.null(drop.rows)) {
       if(! (is.character(drop.rows) | is.integer(drop.rows) | is.numeric(drop.rows) ) )
         stop(" drop.rows, when a vector, must be of class character ,integer, or numeric. ")
  }


 ## check keep.mrks and keep.rows
  if(!is.null(keep.mrks) | !is.null(keep.rows)) {
     if(length(keep.chrm) > 1)
        stop(" When keep.mrks and/or keep.rows are specificed, only a single chromosome can be specified.")
  }

  ## check keep.mrks
  if(!is.null(keep.mrks)) {
       if(!(is.character(keep.mrks) | is.integer(keep.mrks) | is.numeric(keep.mrks)) )
         stop(" keep.mrks, when a vector, must be of class character, integer, or numeric. ")
  }

  ## check keep.rows
  if(!is.null(keep.rows)) {
       if(! (is.character(keep.rows) | is.integer(keep.rows) | is.numeric(keep.rows) ) )
         stop(" keep.rows, when a vector, must be of class character ,integer, or numeric. ")
  }


  ## check that drop and keep rows are both being specified 
  if(!is.null(drop.rows) & !is.null(keep.rows))
        stop(" drop.rows and keep.rows cannot both be specified.")

  ## check that drop and keep mrks are both being specified 
  if(!is.null(drop.mrks) & !is.null(keep.mrks))
        stop(" drop.mrks and keep.mrks cannot both be specified.")

 }  ## end i.checks.subsetdata




i.order <- function(x)
{
 ## internal function
 ## ARgs:
 ##      x  of class dmdatachrm
 ## Output
 ##      object of class "ser_permutation" "list"  

  ## Order rows/lines based on deletion pattern
  d <- x
  class(d) <- "matrix"
  sord <- seriate(d, method="BEA")

  return(sord)

}


i.getorder <- function(x, blocks, sord   )
{
 ## internal function
 ## Args:
 ##    x      an object of class dmdatachrm
 ##    blocks  a named integer vector. Elements of vector give which block
 ##            marker loci belong in. 
 ##    sord   an object of type "ser_permutation" and "list"
 ## Outputs
 ##    list object. Each element a named interger vector. The elements
 ##    of the vector give which marker loci are the same.

 # create results structure
 list.orders <- list(blocks=blocks, orders=list())

 for( jj in unique(blocks) )
  {

     cindx <- which(blocks==jj)
     if(length(cindx)==1)
     {
      list.orders$orders[[jj]] <- 1
      names(list.orders$orders[[jj]]) <- names(blocks)[cindx]
     } else {
     dat <- x[, cindx]
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
     names(list.orders$orders[[jj]]) <- colnames(x)[cindx]
     } ## end if else
    } ## end jj
  return(list.orders)
} ## end internal function


i.dmplot.check <- function(x, chrm, keep.mrks, keep.rows)
{
  ## internal function to do some checks for dmplot
  if (!is.dmdata(x)  ) stop("Object not of class  dmdata")

  if(!is.null(chrm))
  {  ## chrm have been specified 
      ## do some checks
      if (!is.numeric(chrm) & !is.integer(chrm))
         stop(" chrm must be either numeric or an integer vector.")
      if (any(chrm < 1 | chrm > attributes(x)$nchrm))
         stop(" chrm must be an integer vector between 1 and the number of chromosomes.")
  } ## end is.null

  if(!is.null(keep.mrks))
  {
    ## markers have been specified but only makes sense if a single chrm has been specified also
    if(!is.numeric(keep.mrks) & !is.integer(keep.mrks))
         stop(" keep.mrks must be either numeric or an integer vector.")
    if(is.null(chrm) | length(chrm) > 1)
         stop(" A single chromosome must be specified.")
    if(ncol(x[[chrm]]) < max(keep.mrks))
         stop(c(" Marker range cannot extend beyond the number of marker loci on the chromosome which is ", ncol(x[[chrm]])))
  } ## end is.null

  if(!is.null(keep.rows))
  {
    ## markers have been specified but only makes sense if a single chrm has been specified also
    if(!is.numeric(keep.rows) & !is.integer(keep.rows))
         stop(" keep.rows must be either numeric or an integer vector.")
    if(is.null(chrm) | length(chrm) > 1)
         stop(" A single chromosome must be specified.")
    if(nrow(x[[chrm]]) < max(keep.rows))
         stop(c(" Row range cannot extend beyond the number of rows which is ", nrow(x[[chrm]])))
  } ## end is.null


}






   
i.dmapping.tidy <- function(res)
{
     ## internal function for dmapping to tidy up final results

       final <- list()

       ## tmp is needed 
       cl <- with(res,    attr(bestx, "class"))   ## keep
       imp <-  with(res, attr(bestx, "imputed"))  ## keep
       indx <- with(res, which(colnames(bestx)=="cut"))
       res <- within(res, tmp <- bestx[,-indx])
       res <- within(res,  bestx <- tmp)
       res <- within(res, attr(bestx, "imputed") <- imp[,-indx])
       res <- within(res, attr(bestx, "class") <- cl)


       ## only keep certain elements of res
       final[["bestx"]] <- res[["bestx"]]
       final[["bestt"]] <- res[["bestt"]]
       final[["besttvec"]] <- res[["besttvec"]]
       if(!is.null(res[["bestdist"]])) final[["bestdist"]] <- res[["bestdist"]]

       return(final)
} 

####----------------------------------##
##   EXTERNAL FUNCTION DECLARATION  ##
###-----------------------------------##



as.dmdist.list <- function(x)
{
# convert list object into dmdist object. 
# list object must have four elements
# checks
  if(length(x) != 4) stop("List is of incorrect length. Must be of length 4")
  if(!(is.numeric(x[[1]]) & is.numeric(x[[2]]) & is.numeric(x[[3]]) & is.numeric(x[[4]])))  
          stop("Elements of list must be numeric")

   res <- list("nrefblocks"=x[[1]], "nobsblocks"=x[[2]], "mindist"=x[[3]], "minwrong"=x[[4]])
   class(res) <- "dmdist"
   res
}


is.dmdata <- function(x) 
  inherits(x, "dmdata")

is.dmdatachrm <- function(dmat) 
  inherits(dmat, "dmdatachrm")


as.dmdata <- function(x) 
  UseMethod("as.dmdata")

as.dmdatachrm <- function(x) 
  UseMethod("as.dmdatachrm")

as.dmdata.default <- function(x)
{
  if(is.dmdata(x) | is.dmdatachrm(x) | is.list(x))
     x
  else
    {
      stop(" Only objects of type matrix, list, or dmdatachrm can be converted into dmdata objects.")
      NULL
    } 
}


as.dmdatachrm.default <- function(x)
{
  if(is.dmdatachrm(x))
     x
  else
    {
      stop(" Only objects of type matrix can be converted into dmdatachrm objects.")
      NULL
    }

}





is.dmdist <- function(y) 
   inherits(y, "dmdist")

as.dmdist <- function(x) 
   UseMethod("as.dmdist")

as.dmdist.default <- function(x)
{
 if(is.dmdist(x))
    x
 else
    {
    stop(" Object not of list type to be converted into a dmdist object.")
    NULL
    }
}  ## end function



as.dmdatachrm.matrix <- function(x)
{

  # convert matrix object into dmdatachrm object.
  # matrix is given generic marker and line names if col and row names are null.
  # marker loci assumed to belong to the same chromosome
  # the chrmname attribute is added

  # checks
  if(ncol(x) < 2) stop("Matrix must contain more than a single column.")
  if(nrow(x) < 2) stop("Matrix must contain more than a single row.")

  # check if matrix contains probs
  genotypes <- TRUE 
  if(any(x>0 & x<1, na.rm=TRUE)) ## probs
     genotypes <- FALSE


  if(genotypes)
  {
     # convert matrix scores into 0/1 where we assume deletions are least prevalent. 
     mode(x) <- "factor"
     mode(x) <- "numeric" ## trick to give me numbers
     t <- table(x)
     if(length(t)==1) ## only a single value. Will assume this is a nondeletion
        x[which(x== as.numeric(names(t[1])), arr.ind=TRUE)] <- 1
     else {
        if(t[1] > t[2] )
        {
           x[which(x== as.numeric(names(t[1])), arr.ind=TRUE)] <- 0
           x[which(x== as.numeric(names(t[2])), arr.ind=TRUE)] <- 1
        } else {
           x[which(x== as.numeric(names(t[2])), arr.ind=TRUE)] <- 0
           x[which(x== as.numeric(names(t[1])), arr.ind=TRUE)] <- 1
       }  ## end if else
     } ## end if else
  } ## end if genotypes

  ## create generic marker names
  if(is.null(colnames(x))) colnames(x) <- paste("M", 1:ncol(x), sep="")

  ## abbreviate line names
  if(is.null(rownames(x))) rownames(x)  <- paste("L", 1:nrow(x), sep="")

  ## add additional attribute for row,column location of imputed genotypes 
  if(is.null(attributes(x)$imputed))
  {
    if(any(is.na(x)))
          attr(x, "imputed") <- which(is.na(x), arr.ind=TRUE)
     else {
       attr(x,"imputed") <- NA
     }
  }
 ## add genotype attribute
 if(is.null(attributes(x)$genotype))
 {
    attr(x, "genotype") <- FALSE             ## probabilities
    if(genotypes) attr(x, "genotype") <- TRUE ## genotypes
 }

  ## add chrmname attribute in preparation for adding a chromosome name to object
 if(is.null(attributes(x)$chrmname))
    attr(x,"chrmname") <- NA


  ## give class to element of y
  class(x) <- "dmdatachrm"
  return(x)
}




as.dmdata.matrix <- function(x)
{
  # convert matrix object into dmdata object.
  # markers are assumed to belong on a single chromosome 

  y <- vector("list", 1)  # list with a one elements
  names(y) <- "1"  ## chromosome number

  y[[1]] <- as.dmdatachrm(x)  ## assigns structure of class dmdatachrm

  ## adding map element to object
  y[["map"]] <- list("1"= seq(1, ncol(y[[1]])))
 
  ## add attribute for number of chromosomes
  attr(y, "nchrm") <- 1

  # new class of object where object contains matrix data with special attributes
  class(y) <- "dmdata"

  ## making sure attributes are not being lost
  if(is.na(attr(y[[1]], "chrmname")) | is.null(attr(y[[1]], "chrmname")))
       attr(y[[1]], "chrmname") <- "1"

  y
}


as.dmdata.dmdatachrm <- function(x)
{
  # convert dmdatachrm object into dmdata object.
  # markers are assumed to belong on a single chromosome 
  # map object created 

  ynew <- vector("list", 1)  # list with a single element
  names(ynew) <- attributes(x)$chrmname   ## chromosome name

  ## simulate map because none available
  map <- sim.dmmap(0, ncol(x))  ## dummy map
  names(map[[1]]) <- colnames(x)  ## making sure correct marker names are being used
  names(map) <- attributes(x)$chrmname  ## adding the chromosome label to name of map

  ynew[[1]] <- x  ## assigns structure of class dmdatachrm

  ## add attribute for number of chromosomes
  attr(ynew, "nchrm") <- 1

  # new class of object where object contains matrix data with special attributes
  class(ynew) <- "dmdata"

  # add map element to list 
  ynew[["map"]] <- map


  ynew
}


as.dmdata.list <- function(x)
{
  # convert list object with elements of class dmdatachrm into 
  # a dmdata object.

  ## check for map object
  if(!any(names(x)=="map"))
    stop(" This list object does not have a map element.")

  indx <- c(1:length(x))[names(x)!="map"]

  for(ch in indx )
   if(!is.dmdatachrm(x[[ch]])) stop(" The elements of this list must be of class dmdatachrm.")

  chrmname <- rep(NA, length(indx))  ## vector of chromosome names to be extracted from dmdatachrm objects
  for(ch in indx)
   chrmname[which(ch==indx)] <- attributes(x[[ch]])$chrmname

  y <- vector("list", length(x))
  names(y)[indx]  <-  chrmname
  y <- x
  class(y) <- "dmdata"
  attr(y, "nchrm") <- length(indx)
  return(y) 
}


sim.dmdata <- function(numlines=NULL, plines=0.5, map=NULL, numblocks=1, Enumdel=4, 
                         p.missing=0, seed=NULL)
{

  ## Purpose:        to simulate deletion data.
  ## Functionality:  Generates a 2-dimensional array for each chromosome in the marker map.
  ##                 The rows are the plant lines and 
  ##                 the columns are the genetic marker loci. The array contains marker 
  ##                 gentoypes 1/0 where 1 denotes a deletion (absent) and 0 denotes 
  ##                 a non-deletion (present).  
  ##
  ##
  ## Args
  ## numlines   -  number of HIB plants in population
  ## plines     -  prob of a line carrying a deletion
  ## map        -  marker map in list format
  ## numblocks  -  minimum number of blocks per chromosome
  ## Enumdel    -  expected number of deletions for poisson distribution of deletions
  ## p.missing  -  prob of a genotype being missing
  ## seed       -  a numeric value to initialize the pseudo-random number generator
  ##
  ## Returns    -  returns two dmdata objects, one with all the genotypes observed
  ##               and the other with missing genotypes (if p.missing is not null). 

  ##-------------------##
  ##   Checks          ##
  ##-------------------##

  if(is.null(numlines)) stop(" Must specify number of lines.")
  if(!all.equal(numlines, as.integer(numlines))) stop(" numlines must be an integer number.\n")
  if(is.null(map)) cat(" Data will be generated where all marker loci are assumed to be on the same chromosome.\n")
  if(numlines < 1) stop(" Must have an integer number of lines.") 
  if(is.null(map)) {
        cat(" Must either specify a map or specify the number of marker loci. If the number of marker loci \n")
        cat(" is specified, then it is assumed that these loci belong to the same chromosome.\n")
        stop(" sim.dmdata has stopped with errors in the parameter statement.\n")
   }


  if(p.missing <0 | p.missing > 1) stop(" p.missing must be a probability between 0 and 1.\n")
  if(numblocks < 1) stop(" numblocks must be a positive integer.\n")
  if(!all.equal(numblocks, as.integer(numblocks))) stop(" numblocks must be an integer number.\n")
   if(Enumdel < 2) stop(" Enumdel must be a number  greater than 1.\n")

  ## set seed if not null
  if(!is.null(seed)) set.seed(seed)


  dm <- list("missing"= vector("list", length(map)), "nomissing"=vector("list", length(map)))
  names(dm[["missing"]]) <- as.character(1:length(map))
  names(dm[["nomissing"]]) <- as.character(1:length(map))

  ## simulate dummy map if map is integer number
  if(is.integer(map)) map <- sim.dmmap(0, map) 

  for(ch in 1:length(map))
  {  
     # get total number of markers on chromosome
     t.mrk <- length(map[[ch]])

     # divide number of markers along chromosome into numblocks
     t.mrk.block <- round(t.mrk/numblocks)

     if(t.mrk.block <= 1)
     {
        cat(" To few marker loci on chromosome ", names(map)[ch], "to divide into ", numblocks, "blocks.\n")
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))
        stop()
     } 

     dmdat.miss <- NULL    
     dmdat.full <- NULL    
     for(jj in 1: numblocks)
     { 
        if(jj < numblocks){
               indx <- (jj-1)*t.mrk.block +seq(1, t.mrk.block)
         } else {
            indx <- seq( (jj-1)*t.mrk.block+1, t.mrk)
         }
 
        mrks <- names(map[[ch]])[indx]

        ## create 2D data structure
        deldat <- create.structure.SimDeletion(numlines, length(mrks), mrks) 

        ## generate deletions by
        ##  1. sampling lines that contain deletions.
        ##  2. sampling length of deletions.
        ##  3. sampling location of deletions. 
  
        res <- step1.simulation(numlines, plines)

        deldat <- with(res, step23.simulation(rnumlines, Enumdel, numlines, length(mrks), chosen.lines, deldat ))
  
        res <- add.missing.deldat(p.missing, deldat) 

        dmdat.miss <- cbind(dmdat.miss, as.matrix(res[[2]]))
        dmdat.full <- cbind(dmdat.full, as.matrix(res[[1]]))
     } ## end for jj in 1:numblocks 

     dm[["missing"]][[ch]] <- as.dmdatachrm(dmdat.miss)
     attributes(dm[["missing"]][[ch]])$chrmname <- names(map)[ch]
     dm[["nomissing"]][[ch]] <- as.dmdatachrm(dmdat.full)
     attributes(dm[["nomissing"]][[ch]])$chrmname <- names(map)[ch]

  } ## end for
  dm[["missing"]][["map"]] <- map
  dm[["missing"]] <- with(dm, as.dmdata(missing))

  dm[["nomissing"]][["map"]] <- map
  dm[["nomissing"]] <- with(dm, as.dmdata(nomissing))

   return(dm)

}  ## end function  sim.dmdata




print.dmdatachrm <- function(x, ...)
{
  cat(" A summary of this object is: \n\n")
  cat(" The (recoded) deletion data for Chromosome ", attributes(x)$chrmname , " is : \n")
  write.it(x)
  write.mrknames(x)
  write.linenames(x)
} ## end function print.dmdata



print.dmdata <- function(x, ...)
{
  cat(" A summary of this object is: \n\n")
  cat( " Object contains ",attributes(x)$nchrm,"chromosomes. \n")


  indx <- seq(1, length(x))[names(x)!="map"] 
  for(ch in indx)
  {
     cat(" \n\n  CHROMOSOME ", attributes(x[[ch]])$chrmname , "\n")
     cat(" -----------------\n")
     cat(" The (recoded) deletion data are: \n")
     write.it(x[[ch]])
     write.mrknames(x[[ch]])
     write.linenames(x[[ch]])
  }
} ## end function print.dmdata


print.dmdist <- function(x, ...)
{
 cat(" Multivariate Distance Measure: \n\n")
 cat( " No. of blocks         No. of blocks   Total minimum       Number wrongly  \n")
 cat( " from reference        from observed   rearrangement dist  plcaed markers  \n")
 cat( "-------------------------------------------------------------------------  \n")
 cat( "      ",x[[1]], "                   ", x[[2]], "              ", x[[3]], "                ", x[[4]], "\n\n")

} ## end function 


CreateDistMatrix <- function(x=NULL, dmethod="manhattan")
{
## Purpose:  Use Manhattan distance measure on marker loci (columns)
## Args:     matrix or dmdata object
## Return:   distance object for deletion data

   if(is.dmdatachrm(x)){
     class(x) <- "matrix"
   } else {
    if(is.matrix(x))
    {
      x <- as.dmdatachrm(x)
    } else stop(" Object must be of type matrix or dmdatachrm.")
   }
   if(any(is.na(x))) stop("NA values not allowed when calculating distances.")


    return(dist(t(x), method=dmethod))

} ## end function CreateDistMatrix


combinedata <- function(x1=NULL, x2=NULL)
{
 ## Purpose:  cbind two dmdatachrm objects
 ## Return:   new dmdatachrm object 

 if(!is.dmdatachrm(x1)) stop(" First object not of class dmdatachrm.")
 if(!is.dmdatachrm(x1)) stop(" Second object not of class dmdatachrm.")
 if(is.null(x1) & is.null(x2)) return()
 if(is.null(x2))  return(x1)

  ## transpose because there is no cbind.fill in plyr
  res <- data.matrix(plyr:::rbind.fill(as.data.frame(t(x1[[1]][,])), as.data.frame(t(x2[[1]][,]))))
  res <- t(res)
  colnames(res) <- c(colnames(x1[[1]]), colnames(x2[[1]])) 


   PROBLEM - NEED A MAP OBJECT BUT NOt Available with dmdatachrm

}



subsetdata  <- function(x=NULL, keep.chrm=NULL, drop.mrks=NULL, drop.rows=NULL, keep.mrks=NULL, keep.rows=NULL)
{
  ## Purpose: subset dmdata object for markers and/or individuals
  ##          and adjust map accordingly. 
  ## 
  ##          If mrks or individuals are non-NULL, then x is only allowed to
  ##          have a single chromosome
  ##

  ## perform checks
  i.checks.subsetdata(x, keep.chrm, drop.mrks, drop.rows, keep.mrks, keep.rows)
  
  if(length(x)==2) keep.chrm <- 1 ## only a single chromosome 

  ## subsetting on chromosomes
  chrms <- names(x)
  chrms <- chrms[- grep("map", chrms)]  ## removing map label
  keep.chrm <- sort(keep.chrm) ## just in case chrm labels are mis-ordered
  if(!is.null(keep.chrm)) {
     if(is.character(keep.chrm)) {
         indx <- match(keep.chrm, chrms)
         if(any(is.na(indx)))  
             stop(paste(c(" keep.chrm contains the following mismatched chromosome labels ... ", keep.chrm[is.na(indx)]), collapse=" ")) 
     }  ## end if is.character

    if(is.integer(keep.chrm) | is.numeric(keep.chrm)){
      indx <- match(keep.chrm, 1:length(chrms))
      if(any(is.na(indx)))  
             stop(paste(c(" keep.chrm is indexing chromosomes that do not exist ... ", keep.chrm[is.na(indx)]), collapse=" "))
    }  ## end if is.integer | is.numeric

     ## form subsetted structure and adjust map accordingly
     subx <- x[indx]
     attr(subx, "nchrm") <- length(indx)
     indx <- keep.chrm  ## initialize
     if(is.character(keep.chrm)) indx <- match(keep.chrm, names(x[["map"]]))


     if(any(is.na(indx)))  stop(" Problem with matching keep.chrm with map element in dmdata object.")
     
     subx[["map"]] <-  x[["map"]][indx]
     subx <- as.dmdata(subx)

  }  ## end if is.null(keep.chrm)

 ## subsetting on mrks to be removed/dropped
 cindx <- 1:ncol(subx[[1]])
 ## drop markers 
 if(!is.null(drop.mrks)){
      if(!is.character(drop.mrks)){  ## numeric or integer
          if(max(drop.mrks) > ncol(subx[[1]]))  
              stop(paste(c(" drop.mrks is indexing marker columns  that do not exist ... ", max(drop.mrks)),collapse= " "))
          cindx <- cindx[-drop.mrks]
      } else {
       indx <- match(drop.mrks, colnames(subx[[1]]))
       if(any(is.na(indx)))
              stop(paste(c(" drop.mrks contains the following mismatched marker labels ... ", drop.mrks[is.na(indx)]),collapse=" ") )
       cindx <- cindx[-indx]
      } ## end if else
 } ## end if is.null


  ## keep markers
 if(!is.null(keep.mrks)){
    if(!is.character(keep.mrks)){  ## numeric or integer
         if(max(keep.mrks) > ncol(subx[[1]]))
             stop(paste(c(" keep.mrks is indexing marker columns  that do not exist ... ", max(keep.mrks)),collapse= " "))
         cindx <- keep.mrks
     } else {
      indx <- match(keep.mrks, colnames(subx[[1]]))
      if(any(is.na(indx)))
             stop(paste(c(" keep.mrks contains the following mismatched marker labels ... ", keep.mrks[is.na(indx)]),collapse=" ") )
      cindx <- indx
     } ## end if else
 } ## end if is.null

    ## form subsetted structure and adjust map accordingly
 if(!is.null(keep.mrks) | !is.null(drop.mrks))
 {
    map <- subx[["map"]]
    map[[1]]  <- map[[1]][cindx]
    chrmname <- attributes(subx[[1]])$chrmname
    subx[[1]] <- as.dmdatachrm(subx[[1]][, cindx] )
    attr(subx[[1]], "chrmname") <- chrmname ## adding back lost attribute
    subx[["map"]][[1]] <- subx[["map"]][[1]][cindx] 
 } ## end if is.null is.null




 ## subsetting on rows to be removed/dropped
 rindx <- 1:nrow(subx[[1]])
 ## drop rows 
 if(!is.null(drop.rows)){
   if(!is.character(drop.rows)){  ## numeric or integer
       if(max(drop.rows) > nrow(subx[[1]]))
           stop(paste(c(" drop.rows is indexing rows  that do not exist ... ", max(drop.rows)),collapse= " "))
       rindx <- rindx[-drop.rows]
   } else {
    indx <- match(drop.rows, rownames(subx[[1]]))
    if(any(is.na(indx)))
           stop(paste(c(" drop.rows contains the following mismatched marker labels ... ", drop.rows[is.na(indx)]),collapse=" ") )
    rindx <- rindx[-indx]
   } ## end if else
 } ## end if is.null



  ## keep rows
 if(!is.null(keep.rows)){
  if(!is.character(keep.rows)){  ## numeric or integer
       if(max(keep.rows) > nrow(subx[[1]]))
           stop(paste(c(" keep.rows is indexing rows  that do not exist ... ", max(keep.rows)),collapse= " "))
       rindx <- keep.rows
   } else {
    indx <- match(keep.rows, rownames(subx[[1]]))
    if(any(is.na(indx)))
           stop(paste(c(" keep.rows contains the following mismatched marker labels ... ", keep.rows[is.na(indx)]),collapse=" ") )
    rindx <- indx
   } ## end if else
 } ## end if is.null


 ## form subsetted structure and adjust map accordingly
 if(!is.null(keep.rows) | !is.null(drop.rows)){
    chrmname <- attributes(subx[[1]])$chrmname
    subx[[1]] <- as.dmdatachrm(subx[[1]][ rindx, ])
    attr(subx[[1]], "chrmname") <- chrmname ## adding back lost attribute
 }  ## end is.null(drop.mrks) is.null(keep.mrks)


return(subx)


} ## end of function subset.dmdata


 






cleandata <- function(x=NULL,  keep.chrm = NULL, drop.mrks = NULL, drop.rows = NULL,
                               keep.mrks = NULL, keep.rows = NULL, ignoreNA=FALSE,
                               probdel=0.05, verbose=TRUE)
{
  ## Purpose:   cleandata 
  ##            Removing of noninformative  rows (i.e. plants that carry no deletions)  
  ##            If ignoreNA=TRUE, rows are removed without regard to existance of NA's 
  ##            IF ignoreNA=FALSE,  rows are only removed if they do not carry any deletions or NAs. 
  ##  
  ## Works on dmdata but only for a single chromosome.  It will return a cleaned object of the same type
  ## Note:  probdel is  probability at which a deletion is assumed. 
  ##        verbose reports the plants being dropped.
   ## checks
   i.checks.subsetdata(x, keep.chrm, drop.mrks, drop.rows, keep.mrks, keep.rows)


   if(is.null(keep.chrm))
   {
     if(length(x) > 2)
        stop("Object is only allowed to contain a single chromosome worth of data.")
   } 
   if(length(keep.chrm) > 1)
        stop("Object is only allowed to contain a single chromosome worth of data.")

   if(length(x)==2) keep.chrm <- 1

   if(!is.null(keep.chrm))
        x <- subsetdata(x=x, keep.chrm=1, drop.mrks=drop.mrks , 
                       drop.rows=drop.rows , keep.mrks=keep.mrks , 
                       keep.rows=keep.rows )



     indxNA <- which(is.na(x[[1]]), arr.ind = TRUE)

     if(ignoreNA) { ## ignoring NA's 
           indx <- which(rowSums(x[[1]] > probdel , na.rm=TRUE) == 0) 
           if(length(indx) > 0)  x <- subsetdata(x, keep.chrm=1, drop.rows=indx)

     }  else {
            ## taking NA's into account in the count
           indx <- which(rowSums(x[[1]] > probdel, na.rm=TRUE) == 0 )
           indx <- indx[is.na(match(indx, indxNA))] ## only keeps indx's that don't have any NA's in rows
           if(length(indx) > 0) x <- subsetdata(x, keep.chrm=1, drop.rows=indx)
     }  ## end if else     
    if(verbose) {
       cat(" Warning: these rows are being removed ... \n")
       cat(names(indx))
       cat("\n\n\n")

    } 


 #    class(x) <- "dmdatachrm"
 #    if(any(is.na(x)))
 #             attr(x, "imputed") <- which(is.na(x[[ch]]), arr.ind=TRUE)
 #    else {
 #             attr(x, "imputed") <- NA
 #    } ## if else

  return(x)
}  ## cleandata






permutecols <- function(x)
{
  if (!is.dmdata(x)  )
     stop("Object not of class  dmdata")

  # save map and chrmname
  map <- x[["map"]]
  
  # create new structure for results
  res <- vector("list", attributes(x)$nchrm)
  names(res) <- as.character(1:attributes(x)$nchrm)
  # permute cols of each chromosome 
  for(ch in 1:attributes(x)$nchrm)  ## indexed over chromosomes
  {
     cl <- class(x[[ch]])
     chrmname <- attributes(x[[ch]])$chrmname
     ## permute columns
     indx <- sample(1:ncol(x[[ch]]), ncol(x[[ch]]), replace=F)
     a <- x[[ch]][, indx]
     res[[ch]] <- as.dmdatachrm(a)
     attr(res[[ch]], "chrmname") <- chrmname
  }
  # turning res into dmdata object
  class(res) <- "dmdata"
  attr(res, "nchrm") <- attributes(x)$nchrm
  res[["map"]] <- map
  return(res)
}



i.getblocks <- function(x)
{
 ## internal function
 ## Args:  
 ##      x  object of class dmdatachrm
 ## Outputs:
 ##      named integer vector giving blocking indexes

   ## initialization of variables
   marker.groupings <- rep(FALSE, ncol(x)-1)
   res <- rep(NA, ncol(x))
   zeroblockindx <- NULL

   ## Determining separate groupings of markers
   for(ii in 2:ncol(x))
      if (x[,ii-1] %*% x[, ii] > 0)
          marker.groupings[ii-1] <-  TRUE

    if(sum(x[,1])==0) {
       zeroblockindx <- 1
       blockindx <- 1
       res[1] <- zeroblockindx
     } else {
       blockindx <- 1
       res[1] <- blockindx
     }
   for(ii in 2:ncol(x))
   {
      if(sum(x[,ii])==0){
         if(is.null(zeroblockindx)) {
           # zeroblock indx not set yet
           zeroblockindx <- blockindx + 1
           blockindx <- blockindx + 1
          }
         res[ii] <- zeroblockindx
       } else {
        if(!marker.groupings[ii-1])
            blockindx <- blockindx + 1
        res[ii] <- blockindx
      }
   }
   names(res) <- colnames(x)

   return(res)
}  ## end internal function







idmarkerblocks <- function(x)
{
  ## Purpose:  to identify separate blocks of markers. 
  ##           In doing this, I am assuming that the data is 
  #            in the true ordering. Some patterns of deletion data
  ##           can lead to independent blocks of markers being formed.
  ##           These blocks can be in different orders. 
  ## Args:     x  -  a dmdatachrm  or dmdata object  
  ## Returns:  named integer vector of blocks if x is of type dmdatachrm. If 
  ##           x is of type dmdata, then a list with named interger elements is returned. 
  ##
  ## Note:     All columns with no deletions are assigned to be from the same block
  ##           NA's not allowed


  if(any(is.na(x))) stop("Error.  NA's not allowed.")
  if (!is.dmdata(x) & !is.dmdatachrm(x)  )
     stop("Object not of class dmdatachrm or dmdata")

  if( is.dmdatachrm(x))
   {

     res <- i.getblocks(x)

   }  else {
    res <- vector("list", attributes(x)$nchrm)
    names(res) <- as.character(1:attributes(x)$nchrm)

    for(ch in 1:attributes(x)$nchrm)
    {
       res[[ch]] <- i.getblocks(x[[ch]])

    }
  } ## end if else
  return(res)
} ## end function





idmarkersuborder <- function(x)
{

  ## Purpose:  to identify all candidate marker orderings
  ##           for each marker block of loci. 
  ## Args:     object of class dmdatachrm or dmdata
  ##           ord  ordering from seriate function
  ## Returns:  list of lists where each element is a integer vector or possible orderings for each 
  ##           chromosome

  if (!is.dmdata(x) & !is.dmdatachrm(x)  )
     stop("Object not of class dmdatachrm or dmdata")


  if(is.dmdatachrm(x))
  {

    ## Order rows/lines based on deletion pattern
    sord <- i.order(x)

    ## identify any marker blocks
    blocks <- idmarkerblocks(x)

    ## get marker orders within each block. This is 
    ## a list object with the number of elements being 
    ## the number of unique blocks. 
    list.orders <- i.getorder(x, blocks, sord   )
   } else {
   list.orders <- vector("list", attributes(x)$nchrm)
   names(list.orders) <- as.character(1:attributes(x)$nchrm)

   for(ch in 1:attributes(x)$nchrm)
   {
      ## Order rows/lines based on deletion pattern
      sord <- i.order(x[[ch]])

      ## identify any marker blocks
      blocks <- idmarkerblocks(x[[ch]])

      ## get marker orders within each block. This is 
      ## a list object with the number of elements being 
      ## the number of unique blocks. 
      list.orders[[ch]] <- i.getorder(x[[ch]], blocks, sord   )
   }
}

  return(list.orders)

} ## end function idmarkersuborder



#   i.impute <- function(x, uniform=TRUE)
#   {
#    ## internal function
#    ## Args:
#    ##       x object of class dmdatachrm
#     mp <- 0.05
#
#     if(!is.dmdatachrm(x)) stop(" Not of class dmdatachrm.")
#
#
#     if(is.na(attributes(x)$imputed))
#     {  ## not set yet
#     indx <- which(is.na(x), arr.ind=TRUE)
#     attr(x, "imputed") <- indx  ## imputed 0,1 data
#     } else {
#       indx <- attributes(x)$imputed
#     }
#
#     if(uniform)
#        x[indx] <- sample(unique(x)[!is.na(unique(x))], nrow(indx), replace=T)
#
#     if(!uniform) {
#     
#        for(ii in 1:nrow(indx))
#        {
#          ro <- indx[ii,1]; co <- indx[ii,2]
#          if (co > 1 & co < ncol(x) ) {
#             jjr <- jjl <- co
#           
#             while(is.na(x[ro,jjr])  & jjr != ncol(x)) jjr <- jjr+1
#             while(is.na(x[ro,jjl]) &  jjl != 1)           jjl <- jjl- 1
#             if (sum(x[ro, c(jjl, jjr)], na.rm=TRUE) == 0)  x[ro,co] <- 
#                            sample(0:1, 1, prob=c(1-mp, mp))
#             if (sum(x[ro, c(jjl, jjr)], na.rm=TRUE) == 1)  x[ro,co] <- 
#                            sample(0:1, 1, prob=c(0.5, 0.5))
#             if (sum(x[ro, c(jjl, jjr)], na.rm=TRUE) == 2)  x[ro,co] <- 
#                            sample(0:1, 1, prob=c(mp, 1-mp))
#           } ## end if
#
#          if (co==1) {
#             jjr <- co
#             while(is.na(x[ro,jjr])  & jjr != ncol(x)) jjr <- jjr+1
#             if (x[ro, jjr] == 0)  x[ro,co] <- sample(0:1, 1, prob=c(1-mp, mp))
#             if (x[ro, jjr] == 1)  x[ro,co] <- sample(0:1, 1, prob=c(0.5, 0.5))
#           } # end if
#           if (co == ncol(x)) {
#             jjl <- co
#             while(is.na(x[ro,jjl]) & jjl != 1) jjl <- jjl - 1
#             if (x[ro, jjl] == 0)  x[ro,co] <- sample(0:1, 1, prob=c(1-mp, mp))
#             if (x[ro, jjl] == 1)  x[ro,co] <- sample(0:1, 1, prob=c(0.5, 0.5))
#           } # end if
#        } # end for
#     }  ## end if !uniform
#
#     attr(x, "imputed") <- indx
#   #  x <- as.dmdatachrm(x)
#   
#     return(x)
#    } 

i.impute <- function(x)
{
 
#    ## internal function
#    ## Args:
#    ##       x object of class dmdatachrm
#
   if(!is.dmdatachrm(x)) stop(" Not of class dmdatachrm.")

    attr(x, "imputed") <- x[,]  ## sets dimensions for imputed
    if( any( is.na(x)))
    {
       indx <- which(is.na(x), arr.ind=TRUE)
       x[indx] <- 0.05 ## prob of a deletion for NA
    } 
    for(ii in 1:nrow( attr(x, "imputed") ))
    {
      for(jj in 1:ncol( attr(x, "imputed")))
      {
        attr(x, "imputed")[ii,jj] <- sample(0:1,1, prob=c(1-x[ii,jj], x[ii,jj]))

      }
    }

 return(x)
}

impute <- function(x)
{
   if (!is.dmdata(x) & !is.dmdatachrm(x)  )
     stop("Object not of class dmdatachrm or dmdata")

   if(is.dmdatachrm(x))
     res <- i.impute(x)
   else 
   {
      res <- vector("list", attributes(x)$nchrm)
      names(res) <- as.character(1:attributes(x)$nchrm)

      for(ch in 1:attributes(x)$nchrm)
        res[[ch]] <- i.impute(x[[ch]])
      res[["map"]] <- x[["map"]]
      res <- as.dmdata(res)
   } ## end else
   return(res)

}  ## end function 



x=NULL; keep.chrm=NULL; drop.mrks=NULL; drop.rows=NULL
 keep.mrks=NULL; keep.rows=NULL; sumx=NULL; main=NULL
   cx = 0.5; index=NULL


dmplot <- function(x=NULL, keep.chrm=NULL, drop.mrks=NULL, drop.rows=NULL, 
                   keep.mrks=NULL, keep.rows=NULL, sumx=NULL, main=NULL, 
                   cx = 0.5, index=NULL, ...)
{
   ## check inputs
   i.checks.subsetdata(x, keep.chrm, drop.mrks, drop.rows, keep.mrks, keep.rows)
   if(!is.null(main) & !is.character(main)) stop(" main must be a string.") 

   if (is.null(keep.chrm)) keep.chrm <- 1:attributes(x)$nchrm

   ## only a single chromosome
   if(length(x) == 2) keep.chrm <- 1


   ## chrm names
   chrm.names <- names(x)
   if(is.integer(keep.chrm) | is.numeric(keep.chrm)) chrm.names <- names(x)[keep.chrm]
   if(is.character(keep.chrm)) chrm.names <- keep.chrm

   ## subset data
   subx <- x  ## initialize
   if(!is.null(keep.chrm))
      subx <- subsetdata(x, keep.chrm=chrm.names)

   if(!is.null(keep.mrks) | !is.null(keep.rows))  ## only allowed to have a single chromosome
      subx <- subsetdata(subx, keep.chrm=chrm.names, keep.mrks=keep.mrks, keep.rows=keep.rows)


    for(ch in chrm.names)
    {  ## set up graphic window      

       
       ddat <- subx[[ which(ch==chrm.names) ]]

       ## Testing AWG 01/07/2014
       ##library(ggplot2)
       ##library(reshape)

       plant <- markers <- prob.of.del <- NULL ## avoids R CMD check problems  
       longd <- melt(ddat[,])  #long format
       names(longd) <- c("plant","markers","prob.of.del")

       ## changing the factor levels of plant and markers to be in row and column order 
       longd$markers <- factor(longd$markers, levels=colnames(ddat))
       longd$plant <- factor(longd$plant, levels=rownames(ddat))



       p <- ggplot(longd, aes(markers, plant)) + geom_tile(aes(fill=prob.of.del), 
                                                   colour="white") +
                   scale_fill_gradient(low="orange", high="steelblue", limits=c(0,1)) 

      ## ggplot may change the ordering of the rows/columns
      ##gb <- ggplot_build(p)  ## extracts plot settings
      ##xaxislabels <- with(gb, panel$ranges[[1]]$x.labels)
      ##yaxislabels <- with(gb,panel$ranges[[1]]$y.labels )

      if(is.null(sumx)){  ## plot marker labels
        pl <- p +  theme(axis.text.x = element_text(angle=-10, vjust=1)) + xlab("Markers")
      } else {

        if(sumx=="cM") {  ## plot marker locus positions
          cmdist <- signif(subx[["map"]][[ which(ch==chrm.names) ]],  digits=3)
          pl <- p + scale_x_discrete(label=as.numeric(cmdist)) + 
            xlab("cM position") +  theme(axis.text.x = element_text(angle=-10, vjust=1))
        }

        if(sumx=="NAs") {  ## plot proportion of NAs
          pl <- p + scale_x_discrete(label=signif(colSums(is.na(ddat))/nrow(ddat),  digits=3) ) +
            xlab("prop missing obs")  +  theme(axis.text.x = element_text(angle=-10, vjust=1))
         }
         if(sumx=="dels"){ ## plot proportion of deletions (prob > 0.8)
          pl <- p + scale_x_discrete(label=signif(colSums(ddat>0.8 ,na.rm=TRUE)/nrow(ddat),digits=3))  + 
            xlab("prop of lines with del_prob > 80%")  +  theme(axis.text.x = element_text(angle=-10, vjust=1))
         }
    }  ## end if else
   print(pl)
   } ## end for ch


} ## end function dmplot









#===================================================
# Function to collapse like markers into one
# marker column. If a column contains missing values
# i.e. NA then it is not collapsed because don't know
# if NA is 0 or 1.
# Returns a list with the collapsed data and info on which
# columns were collapsed
# dmat: a dmdata object
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
   class(ret) <- "dmdata"
   return(list(data = ret, collapsed = uniq$cols))
}

i.dmappingchrm.check <- function(x, niterates, nwithin, cooling, psampled)
{
  print(c(" psampled = ", psampled))
  ## internal function to check parameters of dmapping
   if (!is.dmdatachrm(x)  )
     stop("Object not of class dmdatachrm or dmdata")
  ## check for integer values
 if (!is.numeric(niterates) & !is.integer(niterates))
    stop(" niterates must be either numeric or an integer vector.")
 if (!is.numeric(nwithin) & !is.integer(nwithin))
    stop(" nwithin  must be either numeric or an integer vector.")


 if(niterates < 1) stop(" niterates is the number of iterations and should be greater than 1.\n")
 if(nwithin < 1) stop(" nwithin is the number of within temperature iterations and should be more than or equal to 1.\n")

 if( cooling < 0 | cooling > 1) stop(" cooling constant must be between 0 and 1.\n")
 if(psampled < 0 | psampled > 1) stop(" psampled  is a probability, hence it must be between 0 and 1\n")

}



i.dmapping.check <- function(x, chrm, niterates, nwithin, cooling, psampled)
{
  ## internal function to check parameters of dmapping
   if (!is.dmdata(x)  )
     stop("Object not of class dmdatachrm or dmdata")

  if(!is.null(chrm))
  {  ## chrm have been specified 
      ## do some checks
      if (!is.numeric(chrm) & !is.integer(chrm))
         stop(" chrm must be either numeric or an integer vector.")
      if (any(chrm < 1 | chrm > attributes(x)$nchrm))
         stop(" chrm must be an integer vector between 1 and the number of chromosomes.")
  } ## end is.null

  ## check for integer values
 if (!is.numeric(niterates) & !is.integer(niterates))
    stop(" niterates must be either numeric or an integer vector.")
 if (!is.numeric(nwithin) & !is.integer(nwithin))
    stop(" nwithin  must be either numeric or an integer vector.")


 if(niterates < 1) stop(" niterates is the number of iterations and should be greater than 1.\n")
 if(nwithin < 1) stop(" nwithin is the number of within temperature iterations and should be more than or equal to 1.\n")

 if( cooling < 0 | cooling > 1) stop(" cooling constant must be between 0 and 1.\n")
 if(psampled < 0 | psampled > 1) stop(" psampled  is a probability, hence it must be between 0 and 1\n")

}


i.initconfig <- function(x, cooling, niterates, nwithin )
{
   ## internal function to initialize configuration
   ##Args 
   ##   x  dmdatachrm
   ##   cooling  temperature cooling value
   ##  niterates the number of iterates across temperature
   res <- list()

   chrm <- attr(x, "chrmname") ## otherwise we loose this attr
   cn <- colnames(x)
   x <- cbind(x, rep(0, nrow(x)))
   colnames(x) <- c(cn , "cut")
 ##   attributes(x)$imputed <- indx  ## adding attr back so that it is not lost
 class(x) <- "dmdatachrm"


   res[["x"]] <- impute(x)  ##  impute any missing genotypes
   attr(res[["x"]], "chrmname") <- chrm
   res[["bestx"]] <- x
   attr(res[["bestx"]], "chrmname") <- chrm
   res[["t"]] <- res[["bestt"]] <- 1000000  # starting tour length.
   res[["bestdist"]] <- list("nrefblocks"=0, "nobsblocks"=0, "mindist"=0, "minwrong"=0)
   res[["besttvec"]] <- rep(NA, niterates*nwithin)

   res[["temp"]] <- 0.5/(cooling**niterates)
   return(res)
}


   i.newconfig <- function(x, method, dmethod)
   {
     ## internal function to create a new configuration
     res <- list()

     ## impute new set of missing marker genotypes
     res[["x"]] <- impute(x)
       
     # find best ordering of new realization and its touring length
     D <- CreateDistMatrix(attr(res[["x"]], "imputed"), dmethod)
     # adjusted D so that cut has 0 distance to every other city
     tmpD <- as.matrix(D)
     indxr <- which(rownames(tmpD)=="cut")
     indxc <- which(colnames(tmpD)=="cut")
     tmpD[indxr,] <- 0
     tmpD[,indxc] <- 0
     D <- as.dist(tmpD)

     tD <- as.TSP(D)  # in TSP format
     ### tD <- insert_dummy(tD, label="cut")  # adding extra dummy city
     solu   <- solve_TSP(tD, method=method)
     res[["t"]] <- attr(solu, "tour_length")
     n.order <- cut_tour(solu, "cut")  # new marker ordering
     res[["order"]] <- c(n.order, ncol(res[["x"]])) # adding "cut" back in at end of ordering
     names(res[["order"]]) <- c(names(n.order), "cut")

     ### change x to new order
     res  <- within(res, x[,] <- x[,order]) # reorder  newly imputed data
     res  <- within(res, colnames(x) <- names(order))
     res  <- within(res, attr(x, "imputed") <- attr(x, "imputed")[, order])
     res  <- within(res, colnames(attr(x, "imputed")) <- names(order))

     ### add in lost attribute
     attr(res[["x"]], "chrmname") <- attr(x, "chrmname")
     return(res) 
}




dmappingchrm <- function( x=NULL, niterates=100, nwithin=100,cooling=0.99 , psampled=0.1,
                            method="concorde", refblockstr=NULL, refmrkord=NULL, dmethod="manhattan",...)
{
 ## core routine for mapping deletions on a single chromosome
i.dmappingchrm.check(x, niterates, nwithin, cooling, psampled)  ## check inputs


 ## create initial configuration to set 
 ## x, bestx, t, bestt, bestdist, tvec, temp
 resold <- i.initconfig(x, cooling, niterates, nwithin )
 temp <- resold[["temp"]]

 counter <- 1
 for(ii in 1:niterates)
 {
     for(jj in 1:nwithin)
     {

       ## create new configuration which returns
       ##  x, order, t
       resnew <- i.newconfig(resold[["x"]], method, dmethod)
       resnew$t
       resnew$order




       # test if we should move to new realization
       rnd <- runif(1,0,1)
       MHprob <- exp( -1*(resnew[["t"]]-resold[["t"]])/temp)

       if (MHprob >= 1 | rnd < MHprob)
       { # accept new realization
           resold[["x"]] <- resnew[["x"]]  # assign newly imputed data to current data
           resold[["t"]] <- resnew[["t"]]
           resold <- within(resold, attr(x, "imputed") <- attr(resnew[["x"]], "imputed"))

       

           #keep track of best realization
           if(resnew[["t"]] < resold[["bestt"]])
           {
                 print(resnew[["t"]])
                 resold[["bestx"]] <- resnew[["x"]]
                 resold[["bestt"]] <- resnew[["t"]]

                 blockstrbest <- idmarkerblocks(resold[["bestx"]])
                 if(!is.null(refblockstr) & !is.null(refmrkord))
                     resold[["bestdist"]] <- RearrangeDist(refblockstr, refmrkord, blockstrbest)
           }
        } ## end if

        if(!is.null(refblockstr) & !is.null(refmrkord))
        {
              obsblockstr <- idmarkerblocks(resold[["oldx"]])
              tmp <- RearrangeDist(refblockstr, refmrkord, obsblockstr) ## calculate distance metric
        }
        resold[["besttvec"]][counter] <- resold[["bestt"]]
        counter <- counter + 1

     }  ## end for inner
     temp <- temp * cooling
  }  ## end for outer

    ##----------------------------##
    ## Tidying of results   ##
    ##----------------------------##
    resold <- i.dmapping.tidy(resold)

    return(resold)

}





dmapping <- function(x=NULL, chrm=NULL, niterates=100, nwithin=100,cooling=0.99 , psampled=0.1,
                            method="concorde", refblockstr=NULL, refmrkord=NULL, dmethod="manhattan", ...)
{
  ## core routine for deletion mapping
    i.dmapping.check(x, chrm, niterates, nwithin, cooling, psampled)  ## check inputs

  ## set up results list
  final <- list( "bestx" = vector("list", attributes(x)$nchrm),
                "bestt" = vector("list", attributes(x)$nchrm),
                "besttvec" = vector("list", attributes(x)$nchrm) )
  if(!(is.null(refblockstr) | is.null(refmrkord) )) final[["bestdist"]] <- vector("list", attributes(x)$nchrm)

  
  if(is.null(chrm)) chrm <- 1:attributes(x)$nchrm

  if(is.null(x[["map"]])) stop(" Object must have map. Map not present...")
  map <- x[["map"]]


  for(ch in chrm)
  {
     res <- dmappingchrm(x[[ch]], niterates, nwithin, cooling, psampled, method, 
                            refblockstr, refmrkord, dmethod, ...)
     final[["bestx"]][[ which(ch==chrm) ]] <- res[["bestx"]]
     final[["bestt"]][[ which(ch==chrm) ]] <- res[["bestt"]]
     final[["besttvec"]][[  which(ch==chrm) ]] <- res[["besttvec"]]
     names(final[["bestx"]])[which(ch==chrm)] <-
     names(final[["bestt"]])[which(ch==chrm)] <- 
     names(final[["besttvec"]])[ which(ch==chrm)] <- attributes(x[[ch]])$chrmname 
     attributes(final[["bestx"]][[ which(ch==chrm) ]])$chrmname  <- attributes(x[[ch]])$chrmname
     ## rearrange the marker loci for chromosome ch to be in the same order as final[["bestx"]][[ch]]
     map <- x[["map"]]
     indx <- match(colnames(final[["bestx"]][[which(ch==chrm)]][,]) ,  names(x[["map"]][[ch]]))
     map[[ch]] <- map[[ch]][indx]  ## reordering chromosome ch


     if(!is.null(refblockstr) & !is.null(refmrkord))
      {
           final[["bestdist"]] <- vector("list", attributes(x)$nchrm)
           final[["bestdist"]][[  attributes(x[[ch]])$chrmname ]] <- res[["bestdist"]] 
      }

  }  ## end for ch
  ## add the map

  final[["bestx"]][["map"]] <- map
  final[["bestx"]] <- as.dmdata(final[["bestx"]])


 return(final)
}  ## end function









##--------------------------------##
## Under Construction             ##
##--------------------------------##

## workign on distance Metric

#d <-   sim.dmdata( numlines = 50, plines = 0.2, nummarkers = 50,
#    Enumdel = 8, p.missing = 0.1,  seed = 1)

#dc <- cleandata(d[["missing"]])


### obtain true blocking structure
#mrks <- colnames(dc)
#cindx <- match(mrks, colnames(d[["nomissing"]]))
#dn <- as.dmdata(d[["nomissing"]][,cindx])
#trueblockstr <- idmarkerblocks(dn)
#truemrkord   <- idmarkersuborder(dn)$orders






RearrangeDist <- function(refblockstr=NULL, refmrkord=NULL, obsblockstr=NULL)
{
## refblockstr    vector object for true block structure (order specific)
## refmrkord      list object of true marker ordering within block structure (order specific)
## obsblockstr vector object of realized block structure from simulated annealing  (order specific)

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
res[["nrefblocks"]] <- length(unique(refblockstr))
res[["nobsblocks"]]  <- length(unique(refblockstr))



##---------------------------------------------##
## Total minumum rearrangement distance        ##
## Total minimum wrongly placed markers        ##
##---------------------------------------------##


## calculate total rearrangement distance
total.dist <- 0  ## initalize total rearrangement distance
total.wrong <- 0 ## initalize total min number of wrongly placed markers
for(ii in unique(obsblockstr))
{
 ## identify which true block best matches the observed block (based on max number of matching markers)
 m.obs <- names(obsblockstr)[which(obsblockstr==ii)] 

 tb <- table(refblockstr[match(m.obs, names(refblockstr))])  ## frequency of markers in true blocks
 best.block <- as.numeric(names(which(max(tb)==tb)))  ## best block (may be more than one if ties

 vdist <- rep(NA, length(best.block))
 vwrong <- rep(NA, length(best.block))  ## wrongly placed markers
 for(jj in best.block)  ## allowing for ties of best block
 {
   ##  only keep matching (intersection of)  markers in m.obs and m
   m <- names(refblockstr)[which(refblockstr==jj)]
   m.obs.match <- m.obs[!is.na(match(m.obs, m))]
   m.match <- m[!is.na(match(m, m.obs))]
   ## recode common markers in m and m.obs but conditional on true marker recoding (mrkord)
   m.obs.recoded <- unlist(refmrkord)[match(m.obs.match, names(unlist(refmrkord)))]  ## using coding of true ordering
   m.recoded <- unlist(refmrkord)[ match(m.match, names(unlist(refmrkord)))]  ## using coding of true ordering
 
   vdist[which(best.block==jj)] <- l1Dist(m.obs.recoded, m.recoded)  ## rearrangement distance
   vwrong[which(best.block==jj)] <- length(m.obs) - length(m.obs.match)  ## number of wrongly placed markers
 }     ## end for jj
 total.dist <- total.dist + min(vdist)
 total.wrong <- total.wrong + min(vwrong)

 
} ## end for ii over unique observed blocks
res[["mindist"]] <- total.dist
res[["minwrong"]] <- total.wrong

res <- as.dmdist(res)

return(res)

} ## end function



##-------------------------##
##  Reading in Data        ##
##-------------------------##


i.rf.checkinputs <- function(datafile, mapfile,   mrks.as.rows, names.pres )
{
  ## internal function to check input parameters for errors
 if(missing(datafile))
    stop(" Data file cannot be missing.")
 if(!file.exists(datafile))
   stop(" Data file does not exist.")
 if(!is.logical(mrks.as.rows))
   stop("mrks.as.rows must be TRUE/FALSE only. If TRUE, the first column contains the names of the markers.")
 if(!is.logical(names.pres))
    stop("names.pres must be TRUE/FALSE only. If TRUE, marker and plant names are present in the file.")

   if(!missing(mapfile))
   {
     if(!file.exists(mapfile))
         stop(" Map  file does not exist.")
   }

}


i.rf.checkdatafile <- function(datafile, names.pres, mrks.as.rows,  sep)
{
  ## internal function to check format of data file for errors
   ## checking for correct number of header labels and fields in records
   cf <- count.fields(datafile, sep=sep)
   if(length(cf)==1)
     stop("File only contains a single line.")

   if (names.pres)
   { ## datafile contains marker and plant names 
     if(length(unique(cf))==2)
        { ## checking to make sure we don't have lines of the same length as header
          if (sum(cf== cf[1])> 1)
          {
             df <- data.frame(row.no=1:length(cf), no.elements=cf)
             cat(" \n Rows have incorrect number of elements \n \n")
             cat(" Number of elements in each row: \n")
             print(df, row.names=FALSE)
             opt <- options(show.error.messages=FALSE)
             on.exit(options(opt))
             stop()
          }  # end sum(cf== cf[1])> 1
        } # end if length(unique(cf))==2

     if (length(unique(cf)) > 2)
     {  ## datafile has incorrect record length 
       df <- data.frame(row.no=1:length(cf), no.elements=cf)
       cat(" \n Rows have incorrect number of elements \n \n")
       cat(" Number of elements in each row: \n")
       print(df, row.names=FALSE)
       opt <- options(show.error.messages=FALSE)
       on.exit(options(opt))
       stop()
     } else { ## datafile has correct record length
      if (cf[1] != cf[2] - 1)
      { ## datafile's header line not of correct length
          if (!mrks.as.rows) stop(" Incorrect number of marker labels ")
          if (mrks.as.rows)  stop(" Incorrect number of plant labels ")
      } # end if cf[1]
    } # end if length(unique(cf))
   }  else {
      ## datafile only contains the data
     if (length(unique(cf))>1)
     { ## datafile has incorrect record length
       df <- data.frame(row=1:length(cf), n.elements=cf)
       cat(" \n Rows have incorrect number of elements \n\n")
       cat(" Number of elements in each row \n")
       cat(c(df, "\n"))
       opt <- options(show.error.messages=FALSE)
       on.exit(options(opt))
       stop()
     } # end if length(unique(cf))
  } ## end if name.pres






}


i.rf.checkmapfile <- function(mapfile, sep)
{
   if(!file.exists(mapfile))
         stop(" Map  file does not exist.")
   ## internal function to check the map file for formatting errors
   ## checking for correct number of header labels and fields in records
   cf <- count.fields(mapfile, sep=sep)
   if(length(cf)==1)
           stop("File only contains a single line.")

   if (length(unique(cf))>1)
   { ## mapfile has incorrect record length
          df <- data.frame(row=1:length(cf), n.elements=cf)
          cat(" \n Rows have incorrect number of elements \n\n")
          cat(" Number of elements in each row \n")
          cat(c(df, "\n"))
          opt <- options(show.error.messages=FALSE)
          on.exit(options(opt))
          stop()
   } # end if length(unique(cf))
}


read.files <- function(datafile, mapfile, na.strings="NA", mrks.as.rows=TRUE, 
                       row.names=1, names.pres=TRUE, genotypes=TRUE,
                       sep="", ... )
{

   i.rf.checkinputs(datafile, mapfile,  mrks.as.rows, names.pres )


   ##------------------------##
   ##  Read data file        ##
   ##------------------------##
   i.rf.checkdatafile(datafile, names.pres, mrks.as.rows,  sep)

   ## read in data file
  if(!names.pres) ## no marker or plant names 
    rt <- read.table(file=datafile, header=FALSE, ...)
  else  ## marker and plant names being supplied
    rt <- read.table(file=datafile, header=TRUE, row.names=row.names, ...)

   # convert rt from data frame to matrix
   rt <- as.matrix(rt)
   if(mrks.as.rows) ## transform
       rt <- t(rt)

  ## Another very important check
  if(genotypes) ## presence/absence
  {
      genos <- as.vector(table(rt))
      if(length(genos) > 2)
      { 
          cat(" Data file has more than two genotype classes: ", names(table(rt)), "\n")
          stop() 
      }
  } else { ## non-deletion probabilities
    if(any(rt<0 | rt>1)) 
    {
      cat(" You have set genotypes=FALSE but datafile contains probabilities not in the range of 0 - 1. \n")
      stop()
    }
    if(any(is.na(rt))) stop("datafile is not allowed to contain NAs when genotypes=FALSE")
  }  ## end if else
#   # clean data
#  rt <- cleandata(rt, datafile, names.pres)
#   
#   # turn into dmdata object 
#   res[["x"]]  <- as.dmdata(rt)

   ##--------------------------##
   ## Read map file if present ##
   ##--------------------------##
   if(!missing(mapfile))
   {
       i.rf.checkmapfile(mapfile, sep)
       mapin <- read.table(file= mapfile, header=TRUE,check.names = FALSE)
       map <- list()
       for (i in names(table(mapin[, 2]))) {
             map[[i]] <- mapin[which(mapin[, 2] == i), 3]
             names(map[[i]]) <- mapin[which(mapin[, 2] == i), 
                 1]
       map[[i]] <- sort(map[[i]])
       }
   } ## end if mising mapfile


   ##------------------------------------##
   ## Form dmdata object from rt matrix  ##
   ##------------------------------------##
   res <- list()
   for(ch in names(map))
   {
     mrks <- names(map[[ch]])
     deldata <- matrix(data=NA, nrow=nrow(rt), ncol= length(mrks))
     indx <- match(mrks, colnames(rt))
     if(any(is.na(indx))){
       cat(" The following marker loci are in the map file but not the data file.\n")
       cat( mrks[is.na(indx)], "\n")
       stop(" Marker names mismatch in map and data files.")
     }
     res[[ch]] <- as.dmdatachrm(as.matrix(rt[, indx]))
     attr(res[[ch]], "chrmname") <- ch
   }
   res[["map"]] <- map
   res <- as.dmdata(res)

   return(res)

}  ## end function





sim.dmmap <- function (len = rep(100, 20), n.mar = 10,  eq.spacing = FALSE) 
{
    if (length(len) != length(n.mar) && length(len) != 1 && length(n.mar) != 
        1) 
        stop("Lengths of vectors len and n.mar do not conform.")
    if (length(len) == 1) 
        len <- rep(len, length(n.mar))
    else if (length(n.mar) == 1) 
        n.mar <- rep(n.mar, length(len))
    n.chr <- length(n.mar)
    map <- vector("list", n.chr)
    names(map) <- as.character(1:n.chr)
    for (i in 1:n.chr) {
       if (!eq.spacing) {
                map[[i]] <- sort(runif(n.mar[i], 0, len[i]))
                map[[i]] <- map[[i]] - min(map[[i]])
       }
       else {
                map[[i]] <- seq(0, len[i], length = n.mar[i] + 
                  1)
                map[[i]] <- map[[i]][-1] - map[[i]][2]/2
            }

        names(map[[i]]) <- paste("D", names(map)[i], "M", 1:n.mar[i], 
            sep = "")

    }
    map
}  ## end function











