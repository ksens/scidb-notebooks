# This program is a rewrite of 1000G.Rmd, in python with a scidb interface library.
#

# does this future proof us for Python-3?
from __future__ import print_function, unicode_literals

import time
import datetime
import numpy as np
import pandas as pd
from scidbpy import connect

# shortcut for accessing Array Functional Language (AFL)
afl = sdb.afl

# convert a delta time in seconds to HH:MM:SS format
def string_of_time(delta_sec):
    return str(datetime.timedelta(seconds=delta_sec))

# connect to scidb
sdb = connect('http://scidb20.vdev.dnanexus.com:8080')


SAMPLE = sdb.wrap_array("KG_SAMPLE")
CHROMOSOME = sdb.wrap_array("KG_CHROMOSOME")
VARIANT = sdb.wrap_array("KG_VARIANT")

# could not figure out how to do this:
#     subset(CHROMOSOME, "chrom='21' or chrom='22'"))
#
# closest I got was:
#     CHROMOSOME[('chrom_id' == 21) | ('chrom_id' == 22)]
#
# This is a workaround:
CHROMOSOME20 = CHROMOSOME[0:2]
VARIANT20 = sdb.merge(VARIANT, CHROMOSOME20)

str(VARIANT20)
VARIANT20.head()

# force evaluation
VARIANT20.eval()

# GENOTYPE <- project(merge(scidb("KG_GENOTYPE_PARSED"),project(VARIANT,'chrom')),
#                    c("allele1","allele2","phased"))

GENOTYPE_PARSED= sdb.wrap_array("KG_GENOTYPE_PARSED")
GENOTYPE_FULL = sdb.merge(GENOTYPE_PARSED, VARIANT20['chrom'])
GENOTYPE = GENOTYPE_FULL[["allele1","allele2","phased"]]
GENOTYPE.att_names

start = time.time()
GENOTYPE.eval()
end = time.time()
print string_of_time(end - start)

# Note: the count function returns an array with the number of non-empty
# cells PER FIELD. This is not the same thing as the number of non-empty cells.
str(SAMPLE)
SAMPLE.count()[0]

str(VARIANT)
VARIANT.count()[0]
VARIANT[0:9,:]

str(GENOTYPE)
start = time.time()
GENOTYPE.count()[0]
end = time.time()
print string_of_time(end - start)

## Transition/transversion ratio

# The transition/transversion ratio (Ti/Tv) is a common quality metric
# for variant call sets. Let's compute Ti/Tv of all the variants with
# respect to the reference genome. This first calculation is on the
# variants only, not the individuals' genotypes, and thus involves
# only a modest amount of data.

# This is work in progress, has not been checked yet
SNP = afl.filter(VARIANT,
                 "(ref=='A' or ref=='G' or ref=='C' or ref=='T') and
                  (alt=='A' or alt=='G' or alt=='C' or alt=='T')")
SNP.count()
SNP.eval()
