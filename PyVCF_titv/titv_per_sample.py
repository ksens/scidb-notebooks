#!/usr/bin/python

# Calculate the ti/tv ratio per person (sample)
#
import time
import vcf
import argparse
import sys
import collections
import numpy

# global variables
ti = 0
tv = 0
num_variants = 0
DBG_INTERVAL = 1000

# Count the number of transversions and inversions per person.
Person = collections.namedtuple('person', ['ti', 'tv'])

# Count ti/tv for all people
#   string -> <ti:int, tv:int>
people = {}

# Convert strings like  0|0, 0|1, 1|0, 1|1 to allele counts
def calc_num_alleles(gt):
    if gt is None or len(gt) != 3:
        return 0
    allele1 = int(gt[0])
    allele2 = int(gt[2])
    if allele1 > 2 or allele2 > 2 or allele1 < 0 or allele2 < 0:
        print("strange value found {}".format(gt))
        exit(1)
        #allele1 = min(1, allele1)
        #allele2 = min(1, allele2)
    result = allele1 + allele2
    # if result > 0:
    #    print(gt, allele1, allele2)
    return result

def is_transition(ref, alt):
    return ((ref == 'A' and alt == 'G') or
            (ref == 'G' and alt == 'A') or
            (ref == 'C' and alt == 'T') or
            (ref == 'T' and alt == 'C'))

def iter_samples(rec, is_trans):
    for s in rec.samples:
        nal = calc_num_alleles(s['GT'])
        if nal == 0:
            continue
        p = people.get(s.sample)
        if p is None:
            p = Person(0, 0)
        if is_trans:
            p_new = Person(p.ti + nal, p.tv)
        else:
            p_new = Person(p.ti, p.tv + nal)
        people[s.sample] = p_new

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
        iter_samples(rec, True)
    else:
        tv += 1
        iter_samples(rec, False)

# This version uses PyVCF functions.
# There are slight calling differences between the DIY version and
# this one.
def account_rec_native(rec):
    global ti, tv, people
    if not rec.is_snp:
        return
    if rec.is_transition:
        iter_samples(rec, True)
        ti += 1
    else:
        iter_samples(rec, False)
        tv += 1

def summarize_people_titv():
    global people
    ti_tv_ratios = []
    for _, p in people.iteritems():
        ratio = float(p.ti) / float(p.tv)
        ti_tv_ratios.append(ratio)
    h = numpy.histogram(ti_tv_ratios)
    print(h)

def print_totals(i):
    global ti, tv
    ratio = float(ti) / float(tv)
    if i is not None:
        print("num records={}  tot={} ti={}  tv={}  ti/tv={}".format(i, num_variants, ti, tv, ratio))
    else:
        print("tot={} ti={}  tv={}  ti/tv={}".format(num_variants, ti, tv, ratio))
    summarize_people_titv()
    sys.stdout.flush()

# Process all records in the file
def process_file(fname, nlimit):
    global num_variants
    print("Processing VCF file {}".format(fname))
    with open(fname, 'r') as f:
        vcf_reader = vcf.Reader(f)
        i = 0
        while True:
            num_variants += 1
            i = i + 1
            if i % DBG_INTERVAL == 0:
                print_totals(i)
            if (nlimit is not None and
                i >= nlimit):
                break
            rec = vcf_reader.next()
            #account_rec_native(rec)
            account_rec(rec)

# parse command line
parser = argparse.ArgumentParser(description='Count ti/tv on a VCF file')
parser.add_argument('--file', dest="filenames", action='append',
                    help='VCF file name')
parser.add_argument('--nlimit', type=int,
                    help='how many records to process per file')
args = parser.parse_args()
if (args.filenames is None or
    len(args.filenames) == 0):
    print("must specify at least one VCF file")
    exit(1)

start_time = time.time()

# process all files
for fname in args.filenames:
    process_file(fname, args.nlimit)

end_time = time.time()
diff = end_time - start_time
print("--- {} seconds ---".format(diff))

# print results
print_totals(None)


