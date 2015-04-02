import vcf

filename = "data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# We have a single-letter biallelic SNP.
# Check if this is a transversion or inversion
def count_snp(ref, alt):
    pass

vcf_reader = vcf.Reader(open(filename, 'r'))
for i in range(0,10):
    rec = vcf_reader.next()
    ref = rec.REF
    alt = rec.ALT
    if (ref is None or
        alt is None):
        continue
    if len(alt) != 1:
        continue
    alt = alt[0]
    if (len(ref) != 1 or len(alt) != 1):
        continue
    # we now have a biallelic SNP
    count_snp(ref, alt)
