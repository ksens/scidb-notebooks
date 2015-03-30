import vcf

filename = "data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"


def read_and_process(rec):
    print rec.REF, rec.ALT

vcf_reader = vcf.Reader(open(filename, 'r'))
for i in range(0,9):
    rec = vcf_reader.next()
    read_and_process(rec)
