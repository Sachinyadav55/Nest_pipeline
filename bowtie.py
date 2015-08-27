#! /project/home/sravishankar9/local/bin/python3.4
import os
import sys
import argparse
import logging
import subprocess
import config as cf
from time import gmtime, strftime


def bowtie(read1,read2):
    logging.info('Aligning fastq files to reference using bowtie')
    base = os.path.splitext(read1)[0]
    outfile = base + '.bam'
    bowtie_args = list()
    if read2 == None:
        bowtie_args = [ cf.bowtie, '-x', cf.index, '-U', read1, '-S',outfile, cf.qual_form,
                cf.mode, '-N', cf.mismatch, '-L', cf.seedlen, '-i', cf.interval, '--n-ceil', cf.ambfunc, '--dpad', cf.dpad, '--gbar', cf.gbar, 
                '--ma', cf.matbonus, '--mp', cf.penrange, '--np', cf.ambpen, '--rdg', cf.gappen, '--rfg', cf.refgappen, '-k', cf.maplim,
                '-D', cf.seedex, '-R', cf.reseed, '-p', cf.threads]
    else :
        bowtie_args = [ cf.bowtie, '-x', cf.index, '-1', read1, '-2', read2, '-S',outfile, cf.qual_form,
                cf.mode, '-N', cf.mismatch, '-L', cf.seedlen, '-i', cf.interval, '--n-ceil', cf.ambfunc, '--dpad', cf.dpad, '--gbar', cf.gbar, 
                '--ma', cf.matbonus, '--mp', cf.penrange, '--np', cf.ambpen, '--rdg', cf.gappen, '--rfg', cf.refgappen, '-k', cf.maplim,
                '-D', cf.seedex, '-R', cf.reseed, '-I', cf.minins, '-X', cf.maxins, cf.orient, '-p', cf.threads]
    print (' '.join(bowtie_args))
    try:
        run_bowtie = subprocess.check_call(' '.join(bowtie_args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        logger.info('Read alignment completed successfully')
    except subprocess.CalledProcessError as ret:
        logger.info('Bowtiefailed with return code: '+ str(ret.returncode))
        sys.exit()
    return(outfile)
 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Bowtie2')
    parser.add_argument('-f', type=str, dest='forward', help='Forward read file path')
    parser.add_argument('-r', type=str, dest='reverse', help='Reverse read file path') 
    args = parser.parse_args()
    logger = logging.getLogger('Nest - Read aligner')
    logger.setLevel(logging.DEBUG)
    stream = logging.StreamHandler()
    stream.setLevel(logging.DEBUG)
    formats = logging.Formatter('%(asctime)s;%(name)s;%(levelname)s;%(message)s')
    stream.setFormatter(formats)
    logger.addHandler(stream)
    ret = bowtie(args.forward, args.reverse)

