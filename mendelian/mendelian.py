#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import vcf
import csv
import time
import numpy as np
import logging
import datetime
import argparse
import subprocess
import configparser
from collections import Counter
from collections import defaultdict
from time import gmtime, strftime

logger = logging.getLogger('Nest - GATK')
logger.setLevel(logging.DEBUG)
stream = logging.StreamHandler()
stream.setLevel(logging.DEBUG)
formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
stream.setFormatter(formats)
logger.addHandler(stream)


def caseselect(vcffile, casepattern, outdir, thresh):
    vcfreader = vcf.Reader(filename=vcffile)
    outfile = vcf.Writer(open('{0}/case_variants.vcf'.format(outdir),'w'),
        vcfreader)
    samples = vcfreader.samples
    groups = defaultdict(list)
    for sample, types in zip(samples, casepattern):
        if types == '+':
            groups['Case'].append(sample)
        elif types == '-':
            groups['Control'].append(sample)
    for lines in vcfreader:
        patient_gt = Counter([lines.genotype(val)['GT'] for val in \
            groups['Case']])
        control_gt = Counter([lines.genotype(val)['GT'] for val in \
            groups['Control']])
        #Filter variants by case=control genotypes
        if control_gt['0/1'] > 0 or control_gt['1/1'] > 0 :
            continue
        if patient_gt['0/1'] + patient_gt['1/1'] < thresh:
            continue
        outfile.write_record(lines)
    return ('{0}/case_variants.vcf'.format(outdir))

def annovar(vcffile, outdir):
    annovar = '{0}/annovar/table_annovar.pl'.format(os.path.dirname(os.path.abspath(__file__)))
    database = 'refGene,ljb26_all,avsnp142,popfreq_all_20150413'
    dbpath = '/data/db/annovar_humandb'
    build = 'hg19'
    operation = 'g,f,f,f'
 
    outfile = os.path.splitext(os.path.basename(vcffile))[0]
    outpath = '{0}/{1}'.format(outdir, outfile)
    ann_args = [annovar, vcffile, dbpath, '-buildver', build, '-out', outpath,
        '-remove', '-protocol', database, '-operation', operation, '-nastring',
        '"."', '-vcfinput']
    try:
        run_anno = subprocess.check_call(' '.join(ann_args), shell=True)
        logger.info('Annovar run completed')
    except subprocess.CalledProcessError as ret:
        logger.info('Annovar failed with return code: {0}'.format(ret.returncode))
        logger.info('Annovar command: {0}'.format(ret.cmd))
        sys.exit()
    return('{0}.hg19_multianno.vcf'.format(outpath))

def dcm(vcffile, outdir, database):
    prots = set([prot[6]+'_'+prot[13] for prot in csv.reader(open(database),delimiter='\t')])
    outfile = os.path.splitext(os.path.basename(vcffile))[0]
    outpath = '{0}/{1}_dcm.vcf'.format(outdir, outfile)
    vcfreader = vcf.Reader(filename=vcffile)
    vcfwriter = vcf.Writer(open(outpath,'w'), vcfreader)
    for lines in vcfreader:
        if None in lines.INFO['AAChange.refGene']:
            lines.INFO['DCM'] = 'False'
        else:
            aac = set([pro.split(':')[0]+'_'+pro.split(':')[-1] for pro in lines.INFO['AAChange.refGene']])
            if aac.intersection(prots):
                lines.INFO['DCM'] = 'True'
            else:
                lines.INFO['DCM'] = 'False'
        vcfwriter.write_record(lines)
    return(outpath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run mendelian filter')
    parser.add_argument('-v', type=str, dest='vcffile', help='Vcf file') 
    parser.add_argument('-o', type=str, dest='outdir', help='Output directory path')
    parser.add_argument('-c', type=str, nargs='+', dest='pattern',
        help='Pattern of case and control in vcf file')
    parser.add_argument('-t', type=int, dest='thresh', help='Case threshold')
    parser.add_argument('-d', type=str, dest='database', help='Database file')
    args = parser.parse_args()
    outfile = caseselect(args.vcffile, args.pattern, args.outdir, args.thresh)
    outfile = annovar(outfile, args.outdir)
    outfile = dcm(outfile, args.outdir, args.database)
 
        
