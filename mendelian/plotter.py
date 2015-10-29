#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import vcf
import time
import numpy as np
import logging
import datetime
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import subprocess
import configparser
import pandas as pd
from collections import Counter
from time import gmtime, strftime

logger = logging.getLogger('Nest - GATK')
logger.setLevel(logging.DEBUG)
stream = logging.StreamHandler()
stream.setLevel(logging.DEBUG)
formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
stream.setFormatter(formats)
logger.addHandler(stream)

def scatterplot(vcffile, outdir):
    vcfreader = vcf.Reader(filename=vcffile)
    samples = vcfreader.samples
    variants = list()
    genotype = {'0/0':'Homozygous_Reference',
                '0/1':'Heterozygous',
                '1/1':'Homozygous_Alternate', 
                '1/2':'Non_Reference_Heterozygous',
                '2/2':'Homozygous_Alternate',
                '0/2':'Heterozygous',
                None :'Homozygous_Reference'}
    lmna = 'chr19:52249680'
    for lines in vcfreader:
        for sam, rank in zip(samples, lines.INFO['Rank']):
            var = dict()
            var['var_id'] = '{0}:{1}'.format(lines.CHROM, lines.POS)
            var['sample'] = sam
            var['refdepth'] = lines.genotype(sam)['AD'][0]
            var['altdepth'] = lines.genotype(sam)['AD'][1]
            var['rank'] = float(rank)
            variants.append(var)
    variants = pd.DataFrame(variants)
    variants.sort(columns='rank', inplace=True)
    sns.set_palette('coolwarm')
    call = variants[(variants.var_id == lmna)]
    xy = [call.refdepth, call.altdepth]
    plt.figure()
    plt.title('Variant ranking')
    sns.lmplot(x='refdepth', y='altdepth', data=variants,
        col='sample', row='rank', hue='rank', palette='coolwarm')
    plt.xlim([0,1000])
    plt.ylim([0,1000])
    plt.savefig('test')
    plt.close()
    plt.figure()
    plt.title('Variant genotypes')
     
    variants.to_csv('test',sep='\t',index_col=0)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run mendelian filter')
    parser.add_argument('-v', type=str, dest='vcffile', help='Vcf file') 
    parser.add_argument('-o', type=str, dest='outdir', help='Output directory path')
    args = parser.parse_args()
    outfile = scatterplot(args.vcffile, args.outdir)
 
        
