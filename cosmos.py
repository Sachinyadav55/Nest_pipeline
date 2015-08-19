import vcf
import os
import csv
import sys
import argparse
import pysam

def cosmic_reader(chrom):
    cosmic_vcf = vcf.Reader(open('CosmicMutantExport.vcf.gz','r'))
    vcf_values = sorted([[records.CHROM,records.POS, records.INFO['GENE'],records.REF,records.ALT[0]] for records in cosmic_vcf.fetch(chrom)])
    return (vcf_values)

def vcf_reader(input_file,chrom):
    vcf_file = vcf.Reader(open(input_file))
    var_records = sorted([[records.CHROM,records.POS,records.REF, records.ALT[0]] for records in vcf_file.fetch('chr'+chrom)])
    return(var_records)

def merge(left, right):
    result = []
    left_idx, right_idx = 0, 0
    while left_idx < len(left) and right_idx < len(right):
        # change the direction of this comparison to change the direction of the sort
        if left[left_idx][1] < right[right_idx][1]:
            left_idx += 1
        elif left[left_idx][1] == right[right_idx][1]:
            if left[left_idx][2] == right[right_idx][3] and left[left_idx][3] == right[right_idx][4] and left[left_idx][0][3:] == right[right_idx][0]:
                left[left_idx].append(right[right_idx][2])
                result.append(left[left_idx])
            left_idx += 1
        else :
            right_idx += 1

    return result

if __name__ == '__main__':
    input_file = sys.argv[1]
    cosmic_csv = csv.writer(open(input_file[:input_file.rfind('.')]+'_onco_mutations.csv','w'),delimiter='\t')
    cosmic_csv.writerow(['CHROM','POSITION','REF','ALT','ONCOGENE'])
    for chrom in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']:
        onco_genes = cosmic_reader(chrom) 
        input_file = sys.argv[1]
        vcf_records = vcf_reader(input_file,chrom)
        cosmic_calls = merge(vcf_records, onco_genes)
        for lines in cosmic_calls:
            cosmic_csv.writerow(lines)
