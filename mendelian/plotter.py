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
                '3/3':'Homozygous_Alternate',
                '4/4':'Homozygous_Alternate',
                '0/2':'Heterozygous',
                '1/3':'Non_Reference_Heterozygous',
                '2/3':'Non_Reference_Heterozygous',
                '3/4':'Non_Reference_Heterozygous',
                '2/4':'Non_Reference_Heterozygous',
                '0/3':'Heterozygous',
                '0/4':'Heterozygous',
                '0/5':'Heterozygous',
                None :'Homozygous_Reference'}
    lmna = 'chr19:52249680'
    for lines in vcfreader:
        for sam, rank in zip(samples, lines.INFO['Rank']):
            var = dict()
            var['var_id'] = '{0}:{1}'.format(lines.CHROM, lines.POS)
            var['sample'] = sam
            var['refdepth'] = lines.genotype(sam)['AD'][0]
            var['altdepth'] = lines.genotype(sam)['AD'][1]
            var['genotype'] = genotype[lines.genotype(sam)['GT']]
            var['rank'] = float(rank)
            if lines.genotype(sam)['GT'] != '0/0' and lines.genotype(sam)['GT'] != None:
                var['type'] = lines.INFO['ExonicFunc.refGene'][0]
            else:
                var['type'] = 'None'
            variants.append(var)
    variants = pd.DataFrame(variants)
    variants.sort(columns='rank', inplace=True)
    sns.set_palette('coolwarm')
    plt.figure()
    plt.title('Variant ranking')
    sns.lmplot(x='refdepth', y='altdepth', data=variants, fit_reg=False,
        col='sample', hue='rank', palette='coolwarm', col_wrap=3,
        size=3, legend=False)
    plt.xlim([0,1000])
    plt.ylim([0,1000])
    plt.savefig('{0}/allele_balance'.format(outdir))
    plt.close()
    plt.figure()
    plt.title('Variant genotypes')
    sns.lmplot(x='refdepth', y='altdepth', data=variants, fit_reg=False,
        col='sample', hue='genotype', col_wrap=3, size=3, legend=True,
        palette=dict(Homozygous_Reference='b', Homozygous_Alternate='r',
                    Heterozygous='g', Non_Reference_Heterozygous='y'))
    plt.xlim([0,1000])
    plt.ylim([0,1000])
    plt.savefig('{0}/genotype'.format(outdir))
    plt.close()
    plt.figure()
    plt.title('Variant Types')
    sns.countplot(x='type', data=variants, hue='sample', palette='Set3')
    plt.xticks(rotation=75)
    plt.gcf().tight_layout()
    plt.savefig('{0}/type'.format(outdir))
    plt.close()
     
    variants.to_csv('{0}/case_variants.tsv'.format(outdir),sep='\t',index_col=0)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run mendelian filter')
    parser.add_argument('-v', type=str, dest='vcffile', help='Vcf file') 
    parser.add_argument('-o', type=str, dest='outdir', help='Output directory path')
    args = parser.parse_args()
    outfile = scatterplot(args.vcffile, args.outdir)
 
        
