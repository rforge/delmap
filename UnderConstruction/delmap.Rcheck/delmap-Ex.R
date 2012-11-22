pkgname <- "delmap"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('delmap')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CleanData")
### * CleanData

flush(stderr()); flush(stdout())

### Name: CleanData
### Title: Clean Data in Readiness for Analysis
### Aliases: CleanData

### ** Examples

# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# Number of marker loci before beign cleaned
nrow(Slist[["missing"]])

# Clean data but where missing genotypes are ignored.
CleanData(Slist[["missing"]])
nrow(CleanData(Slist[["missing"]]))

# Clean data but where colums/rows with missing genotypes are retained
# for later use. 
nrow(CleanData(Slist[["missing"]],ignoreNA=FALSE))



cleanEx()
nameEx("CreateDistMatrix")
### * CreateDistMatrix

flush(stderr()); flush(stdout())

### Name: CreateDistMatrix
### Title: Distance Measure Computation
### Aliases: CreateDistMatrix

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index) for the standard data sets.

# deletion data where 1 is a deletion
dd <- matrix(data=c(0,0,1,1, 
                        0,0,0,1, 
                        1,0,0,0, 
                        0,1,0,0,  
                        1,1,1,0), nrow=5)
# calculate distance matrix for deletion data
CreateDistMatrix(mat=dd)




cleanEx()
nameEx("IdentifyMarkerBlocks")
### * IdentifyMarkerBlocks

flush(stderr()); flush(stdout())

### Name: IdentifyMarkerBlocks
### Title: Identify Marker Blocks
### Aliases: IdentifyMarkerBlocks

### ** Examples

# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# class of object that does not contain any missing data
class(Slist[["nomissing"]])

## Clean data i.e. remove rows/columns without deletion data
Slist[["nomissing"]] <- CleanData(Slist[["nomissing"]])

# print the block and possible marker orderings with each block for
# the simulated deletion data
IdentifyMarkerBlocks(Slist[["nomissing"]])



cleanEx()
nameEx("IdentifyMarkerOrd")
### * IdentifyMarkerOrd

flush(stderr()); flush(stdout())

### Name: IdentifyMarkerOrd
### Title: Identify Marker Orderings
### Aliases: IdentifyMarkerOrd

### ** Examples

# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# class of object that does not contain any missing data
class(Slist[["nomissing"]])

## Clean data i.e. remove rows/columns without deletion data
Slist[["nomissing"]] <- CleanData(Slist[["nomissing"]])

# print the block and possible marker orderings with each block for
# the simulated deletion data
IdentifyMarkerOrd(Slist[["nomissing"]])



cleanEx()
nameEx("ImputeMissingGeno")
### * ImputeMissingGeno

flush(stderr()); flush(stdout())

### Name: ImputeMissingGeno
### Title: Impute Missing Marker Genotypes
### Aliases: ImputeMissingGeno

### ** Examples

# example to be done
# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# class of object that does not contain any missing data
class(Slist[["missing"]])


# impute missing genotypes
ImputeMissingGeno(Slist[["missing"]])






cleanEx()
nameEx("PermuteCols")
### * PermuteCols

flush(stderr()); flush(stdout())

### Name: PermuteCols
### Title: Change Order of Marker Columns
### Aliases: PermuteCols

### ** Examples

# generate deletion data
Slist <-  SimDeletions( numlines = 50, plines = 0.2, nummarkers = 20,
    Enumdel = 4, p.missing = 0.02,  seed = 1)

# print the data, including the marker ordering
print(Slist[["nomissing"]])

# permute the data and print the data, including the new marker ordering
print(PermuteCols(Slist[["nomissing"]]))




cleanEx()
nameEx("ReadData")
### * ReadData

flush(stderr()); flush(stdout())

### Name: ReadData
### Title: Read Deletion Data
### Aliases: ReadData

### ** Examples

## Not run: 
##D # read in deletion data
##D deldat <- ReadData(datafile="./deldata.txt", marker.names=TRUE, line.names=TRUE)
## End(Not run)



cleanEx()
nameEx("SimDeletions")
### * SimDeletions

flush(stderr()); flush(stdout())

### Name: SimDeletions
### Title: Simulation
### Aliases: SimDeletions simulate 'simulate deletion data'

### ** Examples


# To generate a inbred population of 500 plants 
# where genotype data are collected on 20 marker loci,
# the probability of a plant having a genomic deletion is 0.2, 
# the average length of a deletion is 8 marker loci, and 
# the probability of a marker genotype being missing is 0.02,
# use:
S <- SimDeletions( numlines = 500, plines = 0.2, nummarkers = 20,  
    Enumdel = 8, p.missing = 0.02,  seed = NULL)

# print contents of matrix
print(S$nomissing)

# print contents of matrix
print(S$missing)



cleanEx()
nameEx("delmap-package")
### * delmap-package

flush(stderr()); flush(stdout())

### Name: delmap-package
### Title: Deletion Mapping
### Aliases: delmap-package delmap
### Keywords: package

### ** Examples

S <- SimDeletions(numlines = 500, plines = 0.2, nummarkers = 20,
    Enumdel = 8, p.missing = 0.02,  seed = NULL)



cleanEx()
nameEx("print.delmap.data")
### * print.delmap.data

flush(stderr()); flush(stdout())

### Name: print.delmap.data
### Title: Print Deletion Data
### Aliases: print.delmap.data print

### ** Examples

# generate deletion data
S <- SimDeletions( numlines = 500, plines = 0.2, nummarkers = 20,
    Enumdel = 8, p.missing = 0.02,  seed = NULL)

# print contents of data where all genotypes are observed
print(S$nomissing)

#print contents of data where there are missing genotypes
print(S$missing)




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
