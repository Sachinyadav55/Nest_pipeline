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

logger = logging.getLogger('Nest - GATK')
logger.setLevel(logging.DEBUG)
stream = logging.StreamHandler()
stream.setLevel(logging.DEBUG)
formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
stream.setFormatter(formats)
logger.addHandler(stream)

def target(samfile):
    logger.info('Running realigner target creator')
    interval = os.path.splitext(samfile)[0] + '.intervals'
    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    gatk = config['GATK']['gatk']
    maxint = config['GATK']['maxins']
    minreads = config['GATK']['minreads']
    mismatch = config['GATK']['mismatch']
    window = config['GATK']['window']
    threads = config['General']['threads']
    reference = config['General']['reference']
    known = list()
    if config['GATK']['known'] == 'null':
        known = ['-known', config['General']['dbsnp']]
    else:
        for lines in config['GATK']['known'].split(','):
            known.append('-known')
            known.append(lines)
    
    trcr_args = [java, '-jar', gatk, '-T', 'RealignerTargetCreator', '-o', interval,
            '-maxInterval', maxint, '-minReads', minreads, '-mismatch', mismatch,
            '-window', window, '-nt', threads, '-I', samfile, '-R', reference ] + known
    try:
        run_trcr = subprocess.check_call(' '.join(trcr_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Realginer target creator completed')
    except subprocess.CalledProcessError as ret:
        logger.info('GATK failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(interval)


def realigner(samfile, interval):

    logger.info('Running indel realigner')
    bamfile = os.path.splitext(samfile)[0] + '_IR.bam'
    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    gatk = config['GATK']['gatk']
    model = config['GATK']['model']
    lod = config['GATK']['lod']
    entropy = config['GATK']['entropy']
    maxins = config['Bowtie']['maxins']
    maxmoves = config['GATK']['maxmoves']
    maxcon = config['GATK']['maxcon']
    reference = config['General']['reference']
    greedy = config['GATK']['greedy']
    maxreads = config['GATK']['maxreads']
    maxinmem = config['GATK']['maxinmem']
    inre_args = [java, '-jar', gatk, '-T', 'IndelRealigner', '--targetIntervals', interval,
            '-o', bamfile, '-model', model, '-LOD', lod, '-entropy', entropy,
            '--maxConsensuses', maxcon, '-maxIsize', maxins, '-maxPosMove', maxmoves, 
            '-greedy', greedy, '-maxReads', maxreads, '-maxInMemory', maxinmem, '-I', samfile, '-R', reference]
    try:
        run_inre = subprocess.check_call(' '.join(inre_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Indel realignment completed')
    except subprocess.CalledProcessError as ret:
        logger.info('GATK failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(bamfile)

def baserecal(samfile):

    logger.info('Running base recalibrator')
    table = os.path.splitext(samfile)[0] + '.table'
    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    reference = config['General']['reference']
    gatk = config['GATK']['gatk']
    ics =  config['GATK']['ics']
    maxcyc = config['GATK']['maxcyc']
    mcs = config['GATK']['mcs']
    bqsrgpen = config['GATK']['bqsrpen']
    ddq = config['GATK']['ddq']
    idq = config['GATK']['idq']
    threads = config['General']['threads']
    covariates = list()
    for covs in config['GATK']['covariates'].split(','):
        covariates.append('-cov')
        covariates.append(covs)
    knownsites = list()
    if config['GATK']['knownsites'] == 'None':
        knownsites = ['-knownSites',config['General']['dbsnp']]
    else:
        for sites in config['GATK']['knownsites'].split(','):
            knownsites.append('-knownSites')
            knownsites.append(sites)
    bsrc_args = [java, '-jar', gatk, '-T', 'BaseRecalibrator', '-R', reference,
            '-I', samfile, '-o', table, '-ics', ics, '-maxCycle', maxcyc, '-mcs', 
            mcs, '-bqsrBAQGOP', bqsrgpen, '-ddq', ddq, '-idq', idq, '-nct', threads] + covariates + knownsites
    try:
        run_bsrc = subprocess.check_call(' '.join(bsrc_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Base recalibration completed')
    except subprocess.CalledProcessError as ret:
        logger.info('GATK failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(table)

def printreads(samfile, table):

    logger.info('Running print reads')
    bamfile = os.path.splitext(samfile)[0] + '_BQ.bam'
    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    gatk = config['GATK']['gatk']
    reference = config['General']['reference']
    threads = config['General']['threads']
    prre_args = [java, '-jar', gatk, '-T', 'PrintReads', '-R', reference, '-I',
            samfile, '-o', bamfile, '-BQSR', table, '-nct', threads]
    try:
        run_prre = subprocess.check_call(' '.join(prre_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Print reads completed')
    except subprocess.CalledProcessError as ret:
        logger.info('GATK failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(bamfile)

def haplocaller(samfile):

    logger.info('Running haplotype caller')
    vcffile = os.path.splitext(samfile)[0] + '.g.vcf'
    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    gatk = config['GATK']['gatk']
    reference = config['General']['reference']
    dbsnp = config['General']['dbsnp']
    calconf = config['GATK']['calconf']
    emitconf = config['GATK']['emitconf']
    threads = config['General']['threads']
    gtmode = config['GATK']['gtmode']
    hpcl_args = [java, '-jar', gatk, '-T', 'HaplotypeCaller', '-R', reference, '-I', samfile, '--dbsnp', dbsnp,'--emitRefConfidence', 'GVCF', '-stand_call_conf', calconf, '-stand_emit_conf', emitconf, '-o', vcffile, '-nct', threads, '-gt_mode', gtmode]
    try:
        run_hpcl = subprocess.check_call(' '.join(hpcl_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Haplotype caller completed')
    except subprocess.CalledProcessError as ret:
        logger.info('GATk failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(vcffile)

def joint_genotyper(gvcf_list, outdir):

    logger.info('Running joint genotyping')
    vcffile = outdir + '/joint_calls.vcf'
    config = configparser.ConfigParser()
    config.read('config.cfg')
    java = config['General']['java']
    gatk = config['GATK']['gatk']
    reference = config['General']['reference']
    variants = [' --variant '+ val for val in gvcf_list]
    jg_args = [java, '-jar', gatk, '-T', 'GenotypeGVCFs', '-R', reference, '-o', vcffile] + variants
    print (jg_args)
    try:
        run_jg = subprocess.check_call(' '.join(jg_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Joint Genptyper completed')
    except subprocess.CalledProcessError as ret:
        logger.info('GATk failed with return code: {0}'.format(ret.returncode))
        sys.exit()
    return(vcffile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run GATK')
    parser.add_argument('-b', type=str, dest='bamfile', help='Bam file path')
    args = parser.parse_args()
    interval = target(args.bamfile)
    bam = realigner(args.bamfile, interval)
    table = baserecal(args.bamfile)
    bam = printreads(args.bamfile, table)
    vcffile = haplocaller(bam)
