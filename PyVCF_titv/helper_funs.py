def is_transition(ref, alt):
    return ((ref == 'A' and alt == 'G') or
            (ref == 'G' and alt == 'A') or
            (ref == 'C' and alt == 'T') or
            (ref == 'T' and alt == 'C'))

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

