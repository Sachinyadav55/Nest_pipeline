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

def target(samfile):
    logger.info('Running realigner target creator')
    interval = os.path.splitext(samfile)[0] + '.intervals'
    trcr_args = [cf.java, '-jar', cf.gatk, '-T', 'RealignerTargetCreator', '-o', interval,
            '-maxInterval', cf.maxint, '-minReads', cf.minreads, '-mismatch', cf.mismatch,
            '-window', cf.window, '-nt', cf.threads, '-I', samfile, '-R', cf.reference ] + cf.known
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
    inre_args = [cf.java, '-jar', cf.gatk, '-T', 'IndelRealigner', '--targetIntervals', interval,
            '-o', bamfile, '-model', cf.model, '-LOD', cf.lod, '-entropy', cf.entropy,
            '--maxConsensuses', cf.maxcon, '-maxIsize', cf.maxins, '-maxPosMove', cf.maxmoves, 
            '-greedy', cf.greedy, '-maxReads', cf.maxreads, '-maxInMemory', cf.maxinmem, '-I', samfile, '-R', cf.reference]
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
    bsrc_args = [cf.java, '-jar', cf.gatk, '-T', 'BaseRecalibrator', '-R', cf.reference,
            '-I', samfile, '-o', table, '-ics', cf.ics, '-maxCycle', cf.maxcyc, '-mcs', 
            cf.mcs, '-bqsrBAQGOP', cf.bqsrgpen, '-ddq', cf.ddq, '-idq', cf.idq, '-nct', cf.threads] + cf.covariates + cf.knownsites
    print (' '.join(bsrc_args))
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
    prre_args = [cf.java, '-jar', cf.gatk, '-T', 'PrintReads', '-R', cf.reference, '-I',
            samfile, '-o', bamfile, '-BQSR', table, '-nct', cf.threads]
    try:
        run_prre = subprocess.check_call(' '.join(prre_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Print reads completed')
    except subprocess.CalledProcessError as ret:
        logger.info('GATK failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(bamfile)

def haplocaller(samfile):
    logger.info('Running haplotype caller')
    vcffile = os.path.splitext(samfile)[0] + '.vcf'
    hpcl_args = [cf.java, '-jar', cf.gatk, '-T', 'HaplotypeCaller', '-R', cf.reference, '-I', samfile, '--dbsnp', cf.dbsnp, '-stand_call_conf', cf.calconf, '-stand_emit_conf', cf.emitconf, '-o', vcffile, '-nct', cf.threads, '-gt_mode', cf.gtmode]
    try:
        run_hpcl = subprocess.check_call(' '.join(hpcl_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Haplotype caller completed')
    except subprocess.CalledProcessError as ret:
        logger.info('GATk failed with return code: '+ str(ret.returncode))
        sys.exit()
    return(vcffile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run GATK')
    parser.add_argument('-b', type=str, dest='bamfile', help='Bam file path')
    args = parser.parse_args()
    logger = logging.getLogger('Nest - GATK')
    logger.setLevel(logging.DEBUG)
    stream = logging.StreamHandler()
    stream.setLevel(logging.DEBUG)
    formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
    stream.setFormatter(formats)
    logger.addHandler(stream)
    #interval = target(args.bamfile)
    #bam = realigner(args.bamfile, interval)
    table = baserecal(args.bamfile)
    bam = printreads(args.bamfile, table)
    vcffile = haplocaller(bam)
