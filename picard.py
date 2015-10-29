#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import time
import datetime
import argparse
import logging
import subprocess
import configparser
from time import gmtime, strftime

logger = logging.getLogger('Nest - Picard')
logger.setLevel(logging.DEBUG)
stream = logging.StreamHandler()
stream.setLevel(logging.DEBUG)
formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
stream.setFormatter(formats)
logger.addHandler(stream)


def addreadgroup(samfile,base):
    #Start of sorting of bam files

    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    mem = config['General']['mem']
    picard = config['Picard']['picard']
    rgid = config['Picard']['rgid']
    rglb = config['Picard']['rglb']
    rgpl = config['Picard']['rgpl']
    rgpu = config['Picard']['rgpu']
    rgsm = config['Picard']['rgsm']
    rgcn = config['Picard']['rgcn']
    rgds = config['Picard']['rgds']
    rgdt = config['Picard']['rgdt']
    rgpi = config['Picard']['rgpi']
    rgpg = config['Picard']['rgpg']
    rgpm = config['Picard']['rgpm']
    if rgid == 'null':
        rgid = base
        rgsm = base
    logger.info('Adding read group info')
    bamfile = os.path.splitext(samfile)[0] + '_RG.bam'
    addrg_param = [java, '-Xmx'+mem, '-jar', picard, 'AddOrReplaceReadGroups',
        'I='+samfile, 'O='+bamfile, 'SORT_ORDER=coordinate', 'RGID='+rgid, 
        'RGLB='+rglb, 'RGPL='+rgpl, 'RGPU='+rgpu, 'RGSM='+rgsm, 'RGCN='+rgcn,
        'RGDS='+rgds, 'RGDT='+rgdt, 'RGPI='+rgpi, 'RGPG='+rgpg, 'RGPM='+rgpm,
        'CREATE_INDEX=true']
    try:
        run_add = subprocess.check_call(' '.join(addrg_param), shell=True)
        logger.info('Read group information added')
    except subprocess.CalledProcessError as ret:
        logger.info('Picard failed wiht return code: {0}'.format(ret.returncode))
        logger.info('Picard command: {0}'.format(ret.cmd))
        sys.exit()
    return(bamfile)

def markdup(samfile):
    #End of read group addition and start duplicate removal

    logger.info('Marking PCR duplicates')
    bamfile = os.path.splitext(samfile)[0] + '_MD.bam'
    metrics = os.path.splitext(samfile)[0] + '_duplication.metrics'
    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    mem = config['General']['mem']
    picard = config['Picard']['picard']
    dup = config['Picard']['dup']
    dupscore = config['Picard']['dupscore']
    dupregex = config['Picard']['dupregex']
    duppix = config['Picard']['duppix']
    mdup_args = [java, '-Xmx'+mem, '-jar', picard, 'MarkDuplicates',
        'I='+samfile, 'O='+bamfile, 'METRICS_FILE='+metrics,
        'REMOVE_DUPLICATES='+dup, 'ASSUME_SORTED=true',
        'DUPLICATE_SCORING_STRATEGY='+dupscore,
        'READ_NAME_REGEX="'+dupregex+'"',
        'OPTICAL_DUPLICATE_PIXEL_DISTANCE='+duppix,
        'CREATE_INDEX=true']
    try:
        run_mdup = subprocess.check_call(' '.join(mdup_args), shell=True)
        logger.info('PCR duplicates marked/removed')
    except subprocess.CalledProcessError as ret:
        logger.info('Picard failed wiht return code: {0}'.format(ret.returncode))
        logger.info('Picard command: {0}'.format(ret.cmd))
        sys.exit()
    return(bamfile, metrics)

def fixmate(samfile):

    logger.info('Fixing mate pair information')
    bamfile = os.path.splitext(samfile)[0] + '_FM.bam'
    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    mem = config['General']['mem']
    picard = config['Picard']['picard']
    fmate_args = [java, '-Xmx'+mem, '-jar', picard, 'FixMateInformation',
        'I='+samfile, 'O='+bamfile, 'ASSUME_SORTED=true', 
        'ADD_MATE_CIGAR=true','CREATE_INDEX=true']
    try:
        run_fmate = subprocess.check_call(' '.join(fmate_args), shell=True)
        logger.info('Mate information corrected')
    except subprocess.CalledProcessError as ret:
        logger.info('Picard failed wiht return code: {0}'.format(ret.returncode))
        logger.info('Picard command: {0}'.format(ret.cmd))
        sys.exit()
    return(bamfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Picard')
    parser.add_argument('-b', type=str, dest='bamfile', help='Bam file path')
    parser.add_argument('-a', type=str, dest='base', help='Basename')
    args = parser.parse_args()
    #bam = addreadgroup(args.bamfile,args.base)
    bam, metrics = markdup(args.bamfile)
    bam = fixmate(bam)
