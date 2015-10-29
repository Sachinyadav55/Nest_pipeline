#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import vcf
import time
import numpy as np
import logging
import datetime
import argparse
import subprocess
import configparser
from collections import Counter
from time import gmtime, strftime

logger = logging.getLogger('Nest - GATK')
logger.setLevel(logging.DEBUG)
stream = logging.StreamHandler()
stream.setLevel(logging.DEBUG)
formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
stream.setFormatter(formats)
logger.addHandler(stream)

def allelebalance(sample, line):
    call = line.genotype(sample)
    if call.is_het:
        ab = call['AD'][1] / sum(call['AD'])
        return(1- abs(ab - 0.5), ab)
    elif call['GT'] == '0/0':
        try:
            ab = call['AD'][0] / sum(call['AD'])
        except ZeroDivisionError:
            ab = 1
        return(1 - abs(ab -1), ab)
    elif call['GT'] == '1/1':
        try:
            ab = call['AD'][1] / sum(call['AD'])
        except ZeroDivisionError:
            ab = 1
        return(1 - abs(ab - 1), ab)
    elif call['GT'] == None:
        return(1, 1)
    else:
        return(1, 1)

def caseselect(vcffile, outdir):
    vcfreader = vcf.Reader(filename=vcffile)
    outfile = vcf.Writer(open('{0}/case_variants.vcf'.format(outdir),'w'),
        vcfreader)
    samples = vcfreader.samples
    for lines in vcfreader:
        rank = 0
        try:    
            if lines.QUAL >= 30 :
                rank += 1
            if lines.INFO['QD'] >= 2:
                rank += 1
            if lines.INFO['FS'] <= 60:
                rank += 1
            if lines.INFO['MQ'] >= 40:
                rank += 1
            if len(lines.FILTER) > 0:
                rank += 1
            if lines.INFO['MQRankSum'] >= 12.5:
                rank += 1
            if lines.INFO['ReadPosRankSum'] >= 8.0:
                rank += 1
        except (KeyError,TypeError):
            continue
        samples_rank = list()
        allele = list()
        for index, sample in enumerate(samples):
            try:
                ab = allelebalance(sample, lines)
                dp = lines.genotype(sample)['DP']
                sample_rank = (ab[0] * dp)  + rank
            except (TypeError, KeyError, ZeroDivisionError):
                continue
            allele.append('{0:.1f}'.format(ab[1]))
            samples_rank.append('{0:.1f}'.format(sample_rank))
        lines.INFO['Rank'] = ','.join(samples_rank)
        lines.INFO['AB'] = ','.join(allele)
        outfile.write_record(lines)
    return (outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run mendelian filter')
    parser.add_argument('-v', type=str, dest='vcffile', help='Vcf file') 
    parser.add_argument('-o', type=str, dest='outdir', help='Output directory path')
    args = parser.parse_args()
    outfile = caseselect(args.vcffile, args.outdir)
 
        
