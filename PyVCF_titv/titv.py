import vcf
import argparse

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
# There are slight calling differences between the DIY version and
# this one.
def account_rec_native(rec):
    global ti, tv
    if not rec.is_snp:
        return
    if rec.is_transition:
        ti += 1
    else:
        tv += 1

def print_totals(i):
    ratio = float(ti) / float(tv)
    print("num records={}  ti={}  tv={}  ti/tv={}".format(i, ti, tv, ratio))

# parse command line
parser = argparse.ArgumentParser(description='Count ti/tv on a VCF file')
parser.add_argument('--file', dest="filename", help='VCF file name')
args = parser.parse_args()
print("VCF file={}".format(args.filename))

# Process all records in the file
i = 0
with open(args.filename, 'r') as f:
    vcf_reader = vcf.Reader(f)
    while True:
        i = i + 1
        if i % 1000 == 0:
            print_totals(i)
        rec = vcf_reader.next()
        account_rec_native(rec)

print_totals(i)
