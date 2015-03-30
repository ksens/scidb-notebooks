SAMPLE <- scidb("KG_SAMPLE")
CHROMOSOME <- scidb("KG_CHROMOSOME")
VARIANT <- merge(scidb("KG_VARIANT"), subset(CHROMOSOME, "chrom='21' or chrom='22'"))
VARIANT <- scidbeval(VARIANT)

PGENOTYPE <- bind(scidb("KG_GENOTYPE"),"gtf","iif(strlen(gt) != 3,null,gt)")
PGENOTYPE <- bind(PGENOTYPE,"allele1","iif(substr(gtf,0,1) = '.',null,uint8(substr(gtf,0,1)))")
PGENOTYPE <- bind(PGENOTYPE,"allele2","iif(substr(gtf,2,1) = '.',null,uint8(substr(gtf,2,1)))")
PGENOTYPE <- bind(PGENOTYPE,"phased","bool(iif(substr(gtf,1,1) = '|',1,0))")
PGENOTYPE <- project(PGENOTYPE,c("allele1","allele2","phased"))

# This takes about 10 minutes on our SciDB instance
system.time(scidbeval(PGENOTYPE, name="KG_GENOTYPE_PARSED", gc=FALSE))

GENOTYPE <- project(merge(scidb("KG_GENOTYPE_PARSED"),project(VARIANT,'chrom')),
                    c("allele1","allele2","phased"))


# count biallelic SNPs
SNP <- subset(VARIANT, "(ref='A' or ref='G' or ref='C' or ref='T') and
                        (alt='A' or alt='G' or alt='C' or alt='T')")
snps <- count(SNP)
snps


# Hardy-Weinberg equilibrium calculations.
# 
# These four lines take about 10 minutes on our SciDB instance. 
GENOTYPE_SNP <- project(merge(GENOTYPE, SNP, "variant_id"), c("allele1","allele2"))
GT_HOM0_BY_SNP <- aggregate(subset(GENOTYPE_SNP, "allele1=0 and allele2=0"),
                            by="variant_id", FUN="count(*)", eval=TRUE)
GT_HOM1_BY_SNP <- aggregate(subset(GENOTYPE_SNP, "allele1=1 and allele2=1"),
                            by="variant_id", FUN="count(*)", eval=TRUE)
GT_HET_BY_SNP <- aggregate(subset(GENOTYPE_SNP, "(allele1=1 and allele2=0) or
                                                 (allele1=0 and allele2=1)"),
                           by="variant_id", FUN="count(*)", eval=TRUE)

# This is an equivalent calculation. It runs for hours, and does not
# not complete.
GT_CTS <- aggregate(merge(GENOTYPE, SNP, "variant_id"), by=list("variant_id", "allele1", "allele2"), FUN="count(*)")
GT_CTS <- scidbeval(GT_CTS)


