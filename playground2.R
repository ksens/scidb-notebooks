# Ohad,
#
# I meant to show you this earlier re: consolidating the multiple
# aggregation operations into one. It is two different ways of computing
# the single-sample Ti/Tv earlier in the workbook, which is a nice
# example because it runs in <1min. The first version does two separate
# aggregations to count transitions and transversions. The second does
# one aggregation with by="is_transition", to count both quantities in
# one pass. But this seems to take fully twice as much time as the two
# separate passes! Meaning it is actually four times slower walking the
# data.  I did also spend a little bit of time trying to consolidate the
# HWE aggregations, and similarly observed it being super slow, so slow
# I did not wait for it to finish after >1hr or something.  Can you
# verify the same observation? If so, definitely seems like something we
# need to bring up with P4.  Mike

GENOTYPE_HG03209 <- merge(GENOTYPE, SAMPLE$sample_name %==% "HG03209")

tic()
GENOTYPE_HG03209 <- merge(GENOTYPE, SAMPLE$sample_name %==% "HG03209")
GENOTYPE_HG03209_TI <- merge(GENOTYPE_HG03209, SNP$is_transition %==% TRUE, "variant_id")
ti_HG03209 <- aggregate(bind(GENOTYPE_HG03209_TI, "alt_copies", "allele1+allele2"),
                        FUN="sum(alt_copies)")[]
ti_HG03209

GENOTYPE_HG03209_TV <- merge(GENOTYPE_HG03209, SNP$is_transition %==% FALSE, "variant_id")
tv_HG03209 <- aggregate(bind(GENOTYPE_HG03209_TV, "alt_copies", "allele1+allele2"),
                        FUN="sum(alt_copies)")[]
tv_HG03209

ti_HG03209/tv_HG03209
toc()

tic()
aggregate(bind(merge(GENOTYPE_HG03209,project(SNP,"is_transition") ,"variant_id"), "alt_copies", "allele1+allele2"),
                        by="is_transition",
                        FUN="sum(alt_copies)")[]
toc()



# Here was my attempt at consolidating the HWE aggregations, which I
# think took hours. Maybe give it a try before you go to bed or
# something :P
tic()
GT_CTS <- aggregate(merge(GENOTYPE, SNP, "variant_id"),
                    by=list("variant_id", "allele1", "allele2"),
                    FUN="count(allele1)", eval=TRUE)
toc()
