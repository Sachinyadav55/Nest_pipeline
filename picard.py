#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import time
import datetime
import argparse
import logging
import subprocess
import config as cf
from time import gmtime, strftime

def addreadgroup(samfile):
    #Start of sorting of bam files
    logger.info('Adding read group info')
    logger.info(';'.join(['RGID='+cf.rgid, 'RGLB='+cf.rglb, 
            'RGPL='+cf.rgpl, 'RGPU'+cf.rgpu, 'RGSM='+cf.rgsm, 'RGCN='+cf.rgcn, 'RGDS='+cf.rgds, 
            'RGDT='+cf.rgdt, 'RGPI='+cf.rgpi, 'RGPG='+cf.rgpg, 'RGPM='+cf.rgpm]))
    bamfile = os.path.splitext(samfile)[0] + '_RG.bam'
    addrg_param = [cf.java, '-Xmx'+cf.mem, '-jar', cf.picard, 'AddOrReplaceReadGroups', 'I='+samfile,
            'O='+bamfile, 'SORT_ORDER=coordinate', 'RGID='+cf.rgid, 'RGLB='+cf.rglb, 
            'RGPL='+cf.rgpl, 'RGPU='+cf.rgpu, 'RGSM='+cf.rgsm, 'RGCN='+cf.rgcn, 'RGDS='+cf.rgds, 
            'RGDT='+cf.rgdt, 'RGPI='+cf.rgpi, 'RGPG='+cf.rgpg, 'RGPM='+cf.rgpm, 'CREATE_INDEX=true']
    try:
        run_add = subprocess.check_call(' '.join(addrg_param), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Read group information added')
    except subprocess.CalledProcessError as ret:
        logger.info('Picard failed wiht return code: '+ str(ret.returncode))
        sys.exit()
    return(bamfile)

def markdup(samfile):
    #End of read group addition and start duplicate removal
    logger.info('Marking PCR duplicates')
    bamfile = os.path.splitext(samfile)[0] + '_MD.bam'
    metrics = os.path.splitext(samfile)[0] + '_duplication.metrics'
    mdup_args = [cf.java, '-Xmx'+cf.mem, '-jar', cf.picard, 'MarkDuplicates', 'I='+samfile, 'O='+bamfile, 
            'METRICS_FILE='+metrics, 'REMOVE_DUPLICATES='+cf.dup, 'ASSUME_SORTED=true',
            'DUPLICATE_SCORING_STRATEGY='+cf.dupscore, 'READ_NAME_REGEX="'+cf.dupregex+'"',
            'OPTICAL_DUPLICATE_PIXEL_DISTANCE='+cf.duppix,'CREATE_INDEX=true']
    print (' '.join(mdup_args))
    try:
        run_mdup = subprocess.check_call(' '.join(mdup_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('PCR duplicates marked/removed')
    except subprocess.CalledProcessError as ret:
        logger.info('Picard failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(bamfile, metrics)

def fixmate(samfile):
    logger.info('Fixing mate pair information')
    bamfile = os.path.splitext(samfile)[0] + '_FM.bam'
    fmate_args = [cf.java, '-Xmx'+cf.mem, '-jar', cf.picard, 'FixMateInformation',
            'I='+samfile, 'O='+bamfile, 'ASSUME_SORTED=true', 'ADD_MATE_CIGAR=true','CREATE_INDEX=true']
    try:
        run_fmate = subprocess.check_call(' '.join(fmate_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Mate information corrected')
    except subprocess.CalledProcessError as ret:
        logger.info('Picard failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(bamfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Picard')
    parser.add_argument('-b', type=str, dest='bamfile', help='Bam file path')
    args = parser.parse_args()
    logger = logging.getLogger('Nest - Picard')
    logger.setLevel(logging.DEBUG)
    stream = logging.StreamHandler()
    stream.setLevel(logging.DEBUG)
    formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
    stream.setFormatter(formats)
    logger.addHandler(stream)
    bam = addreadgroup(args.bamfile)
    bam, metrics = markdup(bam)
    bam = fixmate(bam)
