#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import time
import logging
import datetime
import argparse
import subprocess
import configparser
from time import gmtime, strftime

logger = logging.getLogger('Nest - FastqFilter')
logger.setLevel(logging.DEBUG)
stream = logging.StreamHandler()
stream.setLevel(logging.DEBUG)
formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
stream.setFormatter(formats)
logger.addHandler(stream)

def run_qc(fwd, rev, outdir, base):
    '''Run FastQC and Trimmomatic, written to perform 
    quality trimming in FASTQ files. In the current implementation of this module
    arguments are directly read from the config. Parameters are set to default.
    And can be changed by editing the config file, or by running the 
    pipeline with the autoconfig script. This module is a part of the NEST 
    pipeline, developed by Shashidhar Ravishankar, at the Vannberg Lab, 
    Georgia Institute of Technolgy.'''
    logger.info('Running Trimmomatic')
    config = configparser.ConfigParser()
    config.read('config.cfg')
    fastqc_path = config['FastQC']['fastqc']
    kmer = config['FastQC']['kmer']
    fastqc_dir = outdir + '/fastqc'
    adapters = config['Trimmomatic']['illuminaclip']
    minlen = config['Trimmomatic']['minlen']
    window = config['Trimmomatic']['window']
    leading = config['Trimmomatic']['leading']
    headcrop = config['Trimmomatic']['headcrop']
    crop = config['Trimmomatic']['crop']
    trim_path = config['Trimmomatic']['trimmomatic']
    trailing = config['Trimmomatic']['trailing']
    java = config['General']['java']
    read1 = fwd
    read2 = rev
    print(read1, read2)
    if not os.path.exists(fastqc_dir):
        os.mkdir(fastqc_dir)
    if read2 == None:        #Single end analysis
        output = outdir + '/' + base + '_trimmed.fq.gz' 
        fastqc_args = [ fastqc, '--extract', '-o', fastqc_dir, '-f', 'fastq', 
            '--quiet', read1]
        trim_args = [java, '-jar', trim_path, 'SE', '-phred33', output, 
                'ILLUMINACLIP:{0}:2:30:10'.format(adapters), read1,
                'LEADING:{0}'.format(leading), 'TRAILING:.{0}'.format(trailing), 
                'SLIDINGWINDOW:{0}'.format(window), 'MINLEN:{0}'.format(minlen)]
        try:
            run_fastqc  =  subprocess.check_call(' '.join(fastqc_args), 
                shell=True)
            logger.info('Fastq quality reports created')
        except subprocess.CalledProcessError as ret:
            logger.info('FastQC fail with return code: \
                {0}'.format(ret.returncode))
            logger.info('FastQC comman: {0}'.format(ret.cmd))
            sys.exit()
        try:
            run_trim = subprocess.check_call(' '.join(trim_args), shell=True)
            logger.info('Fastq trimming completed successfully')
        except subprocess.CalledProcessError as ret:
            logger.info('Trimmomatic failed with return code: \
                {0}'.format(ret.returncode))
            logger.info('Trimmomatic command: {0}'.format(ret.cmd))
            sys.exit()
        read1 = output
    else :      #Paired end analysis
        output_fwd = outdir + '/' + base + '_r1_trimmed.fq.gz' 
        output_rev = outdir + '/' + base + '_r2_trimmed.fq.gz' 
        trash_fwd = outdir + '/' + base + '_r1_unpaired.fq.gz' 
        trash_rev = outdir + '/' + base + '_r2_unpaired.fq.gz' 
        fastqc_args = [fastqc_path, '--extract', '-o', fastqc_dir, '-f', 
            'fastq', read1, read2, '--quiet']
        trim_args = [java, '-jar', trim_path, 'PE', '-phred33', read1,
            read2, output_fwd, trash_fwd, output_rev, trash_rev, 
            'ILLUMINACLIP:{0}:2:30:10'.format(adapters), 
            'LEADING:{0}'.format(leading), 'TRAILING:{0}'.format(trailing), 
            'SLIDINGWINDOW:{0}'.format(window), 'MINLEN:{0}'.format(minlen)]
        try:
            run_fastqc = subprocess.check_call(' '.join(fastqc_args),
                shell=True)
            logger.info('Fastq quality reports created')
        except subprocess.CalledProcessError as ret:
            logger.info('FastQC failed with return code: \
                {0}'.format(ret.returncode))
            logger.info('FastQC command: {0}'.format(ret.cmd))
            sys.exit()
        try:
            run_trim = subprocess.check_call(' '.join(trim_args),
                shell=True)
            logger.info('Fastq trimming completed successfully')
        except subprocess.CalledProcessError as ret:
            logger.info('Trimmomatic failed with return code: \
                {0}'.format(ret.returncode))
            logger.info('Trimmomatic command: {0}'.format(ret.cmd))
            sys.exit()
        read1 = output_fwd
        read2 = output_rev
    return (read1,read2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Bowtie2')
    parser.add_argument('-f', type=str, dest='forward', help='Forward read file path')
    parser.add_argument('-r', type=str, dest='reverse', help='Reverse read file path') 
    parser.add_argument('-o', type=str, dest='outdir', help='Output directory')
    parser.add_argument('-b', type=str, dest='base', help='Basename')
    args = parser.parse_args()
    ret = run_qc(args.forward, args.reverse, args.outdir, args.base)
