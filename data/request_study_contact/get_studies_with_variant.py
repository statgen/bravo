#!/usr/bin/env python2

"""
This script receives a variant_id, like '22-36655744-GAGCCGTGTACCCTGAGACC-G' and returns a list of PIs of studies that contain that variant.
"""

import pysam
import csv

def get_studies_with_variant(variant_id):
    chrom, pos, ref, alt = variant_id.split('-')
    assert chrom in ['X', 'Y'] or 1 <= int(chrom) <= 22
    pos = int(pos)
    assert 0 <= pos <= 1e9
    assert all(base in 'ACTG' for base in ref)
    assert all(base in 'ACTG' for base in alt)

    nwd_ids = []
    vcf_filename = '/net/topmed3/working/hmkang/freeze2/10597.v2/final/topmed_freeze2_10597.chr{}.overlap_removed.svm_pass.genotypes.vcf.gz'.format(chrom)
    with pysam.TabixFile(vcf_filename, parser=pysam.asTuple()) as f:
        header = [line for line in f.header if not line.startswith('##')]
        assert len(header) == 1
        header = header[0].split()
        for variant in f.fetch(chrom, pos-1, pos):
            v_pos = int(variant[header.index('POS')])
            v_ref = variant[header.index('REF')]
            v_alt = variant[header.index('ALT')]

            if v_pos != pos or v_ref != ref or v_alt != alt:
                print('Near match: {!r}'.format([pos,v_pos, ref,v_ref, alt,v_alt]))
            else:
                if len(nwd_ids) > 0:
                    raise variant_id # there can only be one
                for colname, data in zip(header, variant):
                    assert colname.startswith('NWD') == (data in ['0/0', '0/1', '1/1']), (colname, data)
                    if data in ['0/1', '1/1']:
                        nwd_ids.append(colname)

    pis = []
    with open('/net/topmed/incoming/study.reference/study.reference/lookup.table.tab') as f:
        samples_info = csv.DictReader(f, delimiter='\t')
        for sample_info in samples_info:
            if sample_info['NWD_ID'] in nwd_ids:
                pis.append(sample_info['PI_NAME'])

    assert len(pis) == len(nwd_ids)
    return set(pis)

if __name__ == '__main__':
    import sys
    variant_id = sys.argv[1]
    for pi in get_studies_with_variant(variant_id): print(pi)
