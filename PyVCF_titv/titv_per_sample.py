a#!/usr/bin/python

# Calculate the ti/tv ratio per person (sample)
#
import vcf
import argparse
import sys
import collections

# global variables
ti = 0
tv = 0
num_variants = 0

# Count the number of transversions and inversions per person.
Person = collections.namedtuple('person', ['ti', 'tv'])

# Count ti/tv for all people
#   string -> <ti:int, tv:int>
people = {}

# Convert strings like  0|0, 0|1, 1|0, 1|1 to allele counts
def calc_num_alleles(gt):
    if gt is None or len(gt) != 3:
        return 0
    allele1 = 0
    allele2 = 0
    if gt[0] == '1':
        allele1 = 1
    if gt[2] == '1':
        allele2 = 1
    return allele1 + allele2

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
    for s in rec.samples:
        nal = calc_num_alleles(s['GT'])
        p = people.get(s.sample)
        if p is None:
            p = Person(0, 0)
        if rec.is_transition:
            p_new = Person(p.ti + 1, p.tv)
        else:
            p_new = Person(p.tv, p.tv + 1)
        people[s.sample] = p_new

def summarize_people_titv():
    ti_tv_array = []
    for p in people:
        ratio = float(p.ti) / float(p.tv)
        ti_tv_array.append(ratio)

def print_totals(i):
    ratio = float(ti) / float(tv)
    if i is not None:
        print("num records={}  tot={} ti={}  tv={}  ti/tv={}".format(i, num_variants, ti, tv, ratio))
    else:
        print("tot={} ti={}  tv={}  ti/tv={}".format(num_variants, ti, tv, ratio))
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
            if i % 1000 == 0:
                print_totals(i)
            if (nlimit is not None and
                i >= nlimit):
                break
            rec = vcf_reader.next()
            account_rec_native(rec)

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

# process all files
for fname in args.filenames:
    process_file(fname, args.nlimit)

# print results
print_totals(None)
summarize_people_titv()


