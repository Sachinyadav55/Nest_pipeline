import os
import sys
import argparse
import subprocess
import config as cf
from time import gmtime, strftime


def read_aligner(read1,read2,reference,threads,outdir,out_log,err_log, phred,sample,alignment):
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('Bowtie2 started at '+ start_time + '\n')
    print('Bowtie2 started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('Bowtie2 started at '+ start_time + '\n')
    sam = outdir + '/' + sample + '.sam' #os.path.basename(read1)[:os.path.basename(read1).rfind['_']] + '.sam'
    if read2 == None:
        run_bowtie = subprocess.Popen([cf.bowtie ,'-p',threads, alignment, '-x',reference,'-1', read1, '-S', sam, '--'+phred],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        bowtie_status = run_bowtie.communicate()
        bowtie_exitcode = run_bowtie.returncode

    else :
        print([cf.bowtie, '-p',threads,alignment, '-x', reference, '-1',read1, '-2', read2, '-S',sam, '--'+phred])
        run_bowtie = subprocess.Popen([cf.bowtie, '-p',threads,alignment, '-x', reference, '-1',read1, '-2', read2, '-S',sam, '--'+phred], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        bowtie_status = run_bowtie.communicate()
        bowtie_exitcode = run_bowtie.returncode


    out_logger.write(bowtie_status[0] + '\n')
    err_logger.write(bowtie_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if bowtie_exitcode:
        out_logger.write('Bowtie crashed at ' +end_time+' with returncode '+str(bowtie_exitcode)+'. Please check error log for details.\n')
        print('Bowtie crashed at ' +end_time+' with returncode '+str(bowtie_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Bowtie crashed at '+ end_time+'\n')
        err_logger.close()
        return(bowtie_exitcode, sam)

    else:
        out_logger.write('Read alignment step completed at '+ end_time)
        print('Read alignment step completed at '+ end_time)
        err_logger.write('Read alignment step completed at '+ end_time)
        out_logger.close()
        err_logger.close()
        return(bowtie_exitcode, sam)

if __name__ == '__main__':
     parser = argparse.ArgumentParser('This module aligns the short reads to refernce genome using bowtie2')
     parser.add_argument('-f','--foward_read',dest='read1', type=str, help='Foward read fastq file path')
     parser.add_argument('-r','--reverse_read',dest='read2', type=str, help='Reverse read fastq file path')
     parser.add_argument('-o','--outdir',type=str,help='Output directory path')
     parser.add_argument('-q','--reference',type=str, default=cf.hg19, help='Minimum sequence quality threshold')
     parser.add_argument('-n','--threads', type=str, default='2', help='Number of threads allocated')
     parser.add_argument('-p','--phred', type=str, default='phred33', help='Phred score')
     parser.add_argument('-a','--alignment_type', type=str, default='local', choices=['local','global'], help='Mode of alignment')
     parser.add_argument('-s','--sample', type=str, help='Sample name')
     args = parser.parse_args()
     out_log = args.outdir + '/log.out'
     err_log = args.outdir + '/log.err'
     align_mode = {'local':'--very-sensitive-local','global':'--very-sensitive'}
     ret = read_aligner(args.read1, args.read2,args.reference, args.threads,args.outdir, out_log, err_log, args.phred,args.sample, align_mode[args.alignment_type])

