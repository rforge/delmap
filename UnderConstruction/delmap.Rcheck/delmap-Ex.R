pkgname <- "delmap"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('delmap')

assign(".oldSearch", search(), pos = 'CheckExEnv')
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
S <-  SimDeletions( numlines = 500, plines = 0.2, nummarkers = 20,
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
