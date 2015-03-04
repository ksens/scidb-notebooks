## R script to load the sample information data
sample.info <- read.csv(header=TRUE, sep="\t", "20140502_complete_sample_summary.txt")

## We are only interested in the sample and population class
m <- subset(sample.info, select = c("Sample", "Population"))

colnames(m) <- c("sample_name", "Population")

## convert into a scidb array, so we can merge it with KG_SAMPLE
msi <- as.scidb(m)

## intersect the KG_SAMPLE with [msi]
z <- merge(x=msi, y=KG_SAMPLE, by.x="sample_name", by.y="sample_name")

# > count(z)
# [1] 2504

zz <- redimension(z, "<Population:string NULL, sample_name:string> [_sample_id=0:2503,2504,0]")
zz <- redimension(zz, "<Population:string NULL, sample_name:string> [sample_id=0:2503,2504,0]")

## save this array as persistent
#
scidbeval(zz, name="KG_SAMPLE_POPL_CODE", gc=FALSE)
KG_SAMPLE_POPL_CODE <- scidb("KG_SAMPLE_POPL_CODE")

## count(KG_SAMPLE_POPL_CODE)
## [1] 2504

