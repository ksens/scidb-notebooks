## ----- experiment -----

## create a simple 1-dimensional array
##
x <- c(1, 2, 1, 1, 2)
X <- as.scidb(x)
## The equivalent of sum(x == 1), sum(x == 2)
## > count(X %==% 1)
## [1] 3
##> count(X %==% 2)
## [1] 2

## In standard R, we can do a 'table(x)', to get the the contingency table.
## In AFL we can do this:
##   redimension (build(<val:int32>[i=1:5,5,0], iif(i<4,1,2)), <n:uint64 null>[val], count(*) as n)"
## redimension(X, "<n:uint64 null>[val]", "count(*) as n")

## redimension(X, "<n:uint32>[i=0:3]", "count(*) as n")

## redimension(X, "<n:uint32>[i=0:3]", "count(*) as n")
## aggregate(X, FUN=table)

m = scidb("build(<val:int32>[i=0:19,20,0], random()%1000)")
m = bind(m, "x", "random()%3")
m = bind(m, "y", "random()%5")
m = scidbeval(m)

## This command
## z <- redimension(m, "<val:double> [x=0:3,3,0, y=0:5,6,0]")

# This works, but does not sum the collisions
z <- redimension(m, "<val:int32> [x=0:3,4,0, y=0:5,6,0]")

## Both of these do not work
##z <- redimension(m, "<val:int32> [x=0:3,4,0, y=0:5,6,0]", FUN="sum(val)")
##z <- redimension(m, "<val:int32> [x=0:3,4,0, y=0:5,6,0]", FUN="sum")

# This works:
z <- redimension(m, "<s:int64 null> [x=0:3,3,0, y=0:5,6,0]",
                 FUN="sum(val) as s")

# Works, but creates an unbounded array.
# We need an additional command, to get to the sub-array.
AGGR <- aggregate(m, by=list("x", "y"), FUN="sum(val) as v", eval=TRUE)
A <- AGGR[0:3, 0:5][]

> m = scidb("build(<val:double>[i=0:99999,100000,0], random())")
> m = bind(m, "x", "random()%1000")
> m = bind(m, "y", "random()%1000")
> m = redimension(m, "<val:double> [x=0:999,1000,0, y=0:999,1000,0]")
> m = scidbeval(m, temp=TRUE)
> count(m)

## ----- experiment -----
# convert the string representation of the phasing into an integer.
#  source string     int
#  0|0               1
#  1|0, 1|0          2
#  1|1               3
#  else              0
X_GENOTYPE <- bind(KG_GENOTYPE, "gti",
                   "int8(iif(gt='0|0', 1, iif(gt='1|1', 3, iif(gt='0|1', 2, iif(gt='1|0', 2, 0)))))")
X_GENOTYPE <- project(X_GENOTYPE, "gti")
system.time(X_GENOTYPE <- scidbeval(X_GENOTYPE), gc=FALSE)
   user  system elapsed
  0.154   0.036 145.408

# These take forever
#   system.time(AGGR <- aggregate(X_GENOTYPE, by=list("variant_id", "gti"), FUN="count(*)", eval=TRUE))
#   system.time(AGGR <- aggregate(X_GENOTYPE, by=list("variant_id", "gti"), FUN="sum", eval=TRUE))
GENOTYPE_SNP <- project(merge(GENOTYPE, SNP, "variant_id"), c("allele1","allele2"))

# Drop the extra chrom_id dimensions
GENOTYPE_SNP_MIN <- redimension(GENOTYPE_SNP,
"<allele1:uint8 NULL DEFAULT null,allele2:uint8 NULL DEFAULT null> [variant_id=0:2229999,10000,0,sample_id=0:2599,100,0]")

# m = redimension(m, "<val:double> [x=0:999,1000,0, y=0:999,1000,0]")
