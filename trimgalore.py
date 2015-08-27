#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import time
import logging
import datetime
import argparse
import subprocess
import config as cf
from time import gmtime, strftime

def run_qc():
    '''Trim galore is a wrapper around fastqc and cutadapt, written to perform 
    quality trimming in FASTQ files. In the current implementation of this module
    arguments are directly read from the config. Parameters are set to default.
    And can be changed by editing the config file, or by running the pipeline with
    the autoconfig script. This module is a part of the NEST pipeline, developed
    by Shashidhar Ravishankar, at the Vannberg Lab, Georgia Institute of Technolgy.'''
    logger.info('Running TrimGalore')
    read1 = cf.read1
    read2 = cf.read2
    if cf.read2 == None:        #Single end analysis
        fastqc_args = [ '--fastqc_args', '"', '-o', cf.outdir, '-t', cf.threads, '-f', 'fastq', 
                '-j', cf.java, '-k', cf.fqc_kmer,'"']
        adapter = open(cf.adapters).readline().strip('\n')
        trim_args = [cf.trim,  cf.qual_form, '-q', cf.qual, '--fastqc']+ fastqc_args + ['-a',adapter, 
                '-e', cf.errorrate, '--stringency', cf.stringency, '--length', cf.length,
                '-o', cf.outdir, '--clip_R1', cf.fiveclip, '--three_prime_clip_R1', cf.threeclip, cf.read1]
        try:
            run_trim = subprocess.check_call(' '.join(trim_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            logger.info('Fastq trimming completed successfully')
        except subprocess.CalledProcessError as ret:
            logger.info('TrimGalore failed with return code: ' + str(ret.returncode))
            sys.exit()
        base = os.path.splitext(cf.read1)[0]
        cf.read1 = base + '_trimmed.fq'

    else :      #Paired end analysis
        fastqc_args = [ '--fastqc_args', '"', '-o', cf.outdir, '-t', cf.threads, '-f', 'fastq', 
                '-j', cf.java, '-k', cf.fqc_kmer,'"']
        adapter = open(cf.adapters).readline().strip('\n')
        trim_args = [cf.trim,  cf.qual_form, '-q', cf.qual, '--fastqc']+ fastqc_args + ['-a',adapter, 
                '-a2', adapter, '-e', cf.errorrate, '--stringency', cf.stringency, '--length', cf.length,
                '-o', cf.outdir, '--clip_R1', cf.fiveclip, '--clip_R2', cf.fiveclip, '--three_prime_clip_R1', 
                cf.threeclip, '--three_prime_clip_R2', cf.threeclip, '--paired', cf.read1, cf.read2]
        try:
            run_trim = subprocess.check_call(' '.join(trim_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            logger.info('Fastq trimming completed successfully')
        except subprocess.CalledProcessError as ret:
            logger.info('TrimGalore failed with return code: '+ str(ret.returncode))
            sys.exit()
        base = os.path.splitext(cf.read1)[0]
        read1  = base + '_trimmed.fq'
        base = os.path.splitext(cf.read2)[0]
        read2 = base + '_trimmed.fq'
    return (read1,read2)

if __name__ == '__main__':
     logger = logging.getLogger('Nest - FastqFilter')
     logger.setLevel(logging.DEBUG)
     stream = logging.StreamHandler()
     stream.setLevel(logging.DEBUG)
     formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
     stream.setFormatter(formats)
     logger.addHandler(stream)
     ret = run_qc()
