import vcf
import pprint

filename = "data/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

ti = 0
tv = 0

def is_transition(ref, alt):
    if ((ref == 'A' and alt == 'G') or
        (ref == 'G' and alt == 'A') or
        (ref == 'C' and alt == 'T') or
        (ref == 'T' and alt == 'C')):
        return True
    else:
        return False

def account_rec(rec):
    global ti, tv
    ref = rec.REF
    alt = rec.ALT
    if (ref is None or
        alt is None):
        return
    if not type(alt) is list:
        return
    if len(alt) != 1:
        return
    alt = alt[0]
    if not type(alt) is vcf.model._Substitution:
        return
    if (len(ref) != 1 or len(alt) != 1):
        return

    # we now have a biallelic SNP
    if is_transition(ref, alt):
        ti += 1
    else:
        tv += 1

# This version uses PyVCF functions.
# There are slight differences, nothing major.
def account_rec_native(rec):
    global ti, tv
    if not rec.is_snp:
        return
    if rec.is_transition:
        ti += 1
    else:
        tv += 1

with open(filename, 'r') as f:
    vcf_reader = vcf.Reader(f)
    #for i in range(0,10000):
    while True:
        rec = vcf_reader.next()
        #account_rec(rec)
        account_rec_native(rec)

print("ti={}  tv={}".format(ti, tv))
