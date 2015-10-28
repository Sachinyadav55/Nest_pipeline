#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import argparse
import logging
import subprocess
import configparser
from time import gmtime, strftime

logger = logging.getLogger('Nest - Read aligner')
logger.setLevel(logging.DEBUG)
stream = logging.StreamHandler()
stream.setLevel(logging.DEBUG)
formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
stream.setFormatter(formats)
logger.addHandler(stream)

def bowtie(read1,read2,outdir,base):

    logger.info('Aligning fastq files to reference using bowtie')
    outfile = outdir + '/' + base + '.bam'
    bowtie_args = list()
    config = configparser.ConfigParser()
    config.read('config.cfg')
    bowtie = config['Bowtie']['bowtie']
    index = config['Bowtie']['index']
    qual_form = '--phred33'
    ambfunc = config['Bowtie']['ambfunc']
    dpad = config['Bowtie']['dpad']
    gbar = config['Bowtie']['gbar']
    matbonus = config['Bowtie']['matbonus']
    penrange = config['Bowtie']['penrange']
    ambpen = config['Bowtie']['ambpen']
    gappen = config['Bowtie']['gappen']
    refgappen = config['Bowtie']['refgappen']
    minins = config['Bowtie']['minins']
    maxins = config['Bowtie']['maxins']
    orient = config['Bowtie']['orient']
    threads = config['General']['threads']
    print(outfile)
    try:
        preset = config['Bowtie']['preset']
        bowtie_args = [ bowtie, '-x', index, '-S',outfile, qual_form,
            preset, '--n-ceil', ambfunc, '--dpad', dpad, '--gbar', gbar, 
            '--ma', matbonus, '--mp', penrange, '--np', ambpen, '--rdg', 
            gappen, '--rfg', refgappen, '-p', threads]
    except KeyError:
        seedex = config['Bowtie']['seedex']
        reseed = config['Bowtie']['reseed']
        seedlen = config['Bowtie']['seedlen']
        mode = config['Bowtie']['mode']
        interval = config['Bowtie']['interval']
        mismatch = config['Bowtie']['mismatch']
        bowtie_args = [ bowtie, '-x', index, '-S',outfile, qual_form,
            mode, '-N', mismatch, '-L', seedlen, '-i', interval, '--n-ceil', 
            ambfunc, '--dpad', dpad, '--gbar', gbar, '--ma', matbonus,
            '--mp', penrange, '--np', ambpen, '--rdg', gappen, '--rfg', 
            refgappen, '-D', seedex, '-R', reseed, '-p', threads]
    try:
        maplim = config['Bowtie']['maplim']
        bowtie_args += ['-k', maplim]
    except KeyError:
        maplim = None
    if read2 == None:
        bowtie_args += ['-U', read1, '--quiet']
    else :
        bowtie_args += ['-1', read1, '-2', read2, '--quiet']
    try:
        run_bowtie = subprocess.check_call(' '.join(bowtie_args), shell=True)
        logger.info('Read alignment completed successfully')
    except subprocess.CalledProcessError as ret:
        logger.info('Bowtie failed with return code: {0}'.format(ret.returncode))
        logger.info('Bowtie command: {0}'.format(ret.cmd))
        sys.exit()
    return(outfile)
 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Bowtie2')
    parser.add_argument('-f', type=str, dest='forward', help='Forward read file path')
    parser.add_argument('-r', type=str, dest='reverse', help='Reverse read file path') 
    parser.add_argument('-o', type=str, dest='outdir', help='Output directory')
    parser.add_argument('-b', type=str, dest='base', help='Basename')
    args = parser.parse_args()
    ret = bowtie(args.forward, args.reverse, args.outdir, args.base)
