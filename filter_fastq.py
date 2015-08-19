import os
import sys
import time
import datetime
import argparse
import subprocess
import config as cf
from time import gmtime, strftime

def run_qc(read1, read2, outdir,threads, out_log, err_log, phred, qual, sample):
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('Fastqc started at '+ start_time + '\n')
    print('Fastqc started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('Fastqc started at '+ start_time + '\n')
    if read2 == None:
        run_fastqc = subprocess.Popen([cf.fastqc,'-t',threads, '-f', 'fastq','-o', outdir, read1],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        fastqc_status = run_fastqc.communicate()
        fastqc_exitcode = run_fastqc.returncode

    else :
        run_fastqc = subprocess.Popen([cf.fastqc,'-t',threads, '-f','fastq', '-o', outdir, read1, read2], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        fastqc_status = run_fastqc.communicate()
        fastqc_exitcode = run_fastqc.returncode
    

    out_logger.write(fastqc_status[0] + '\n')
    err_logger.write(fastqc_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if fastqc_exitcode:
        out_logger.write('Fastqc crashed at ' +end_time+' with returncode '+str(fastqc_exitcode)+'. Please check error log for details.\n')
        print('Fastqc crashed at ' +end_time+' with returncode '+str(fastqc_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Fastqc crashed at '+ end_time+'\n')
        err_logger.close()
        return(fastqc_exitcode, read1, read2)

    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger.write('Sickle started at '+ start_time + '\n')
    print('Sickle started at '+ start_time)
    err_logger.write('Sickle started at '+ start_time + '\n')     
    if read2 == None:
        trim1 = outdir +'/'+sample[0]+'_trimmed.fastq'
        run_sickle = subprocess.Popen([cf.sickle,'se','-f', read1, '-t', phred, '-q', qual,'-o', trim1], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        sickle_status = run_sickle.communicate()
        sickle_exitcode = run_sickle.returncode
   
    else:
        trim1 = outdir + '/' + sample[0] + '_trimmed.fastq' #os.path.basename(read1)[:os.path.basename(read1).rfind('.')]
        trim2 = outdir + '/' + sample[1] + '_trimmed.fastq'
        unp = outdir + '/' + sample[0][:-1] +'unpaired.fastq'
        run_sickle = subprocess.Popen([cf.sickle,'pe','-f',read1,'-r',read2,'-o',trim1,'-p',trim2,'-s',unp,'-t',phred,'-q',qual], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        sickle_status = run_sickle.communicate()
        sickle_exitcode = run_sickle.returncode

    
    out_logger.write(sickle_status[0] + '\n')
    err_logger.write(sickle_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if sickle_exitcode:
        out_logger.write('Sickle crashed at ' +end_time+' with returncode '+str(fastqc_exitcode)+'. Please check error log for details.\n')
        print('Sickle crashed at ' +end_time+' with returncode '+str(fastqc_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Sikcle crashed at '+ end_time +'\n')
        err_logger.close()
        return(sickle_exitcode, trim1, trim2)

    else:
        out_logger.write('Fastq quality check and trimming completed at '+end_time)
        out_logger.close()
        err_logger.write('Fastq quality check and trimming completed at '+end_time)
        err_logger.close()
        print('Fastq quality check and trimming completed at '+end_time)
        return(sickle_exitcode, trim1, trim2)

if __name__ == '__main__':
     parser = argparse.ArgumentParser('This module takes in fastq files as input and perform quality check on them')
     parser.add_argument('-f','--foward_read',dest='read1', type=str, help='Foward read fastq file path')
     parser.add_argument('-r','--reverse_read',dest='read2', type=str, help='Reverse read fastq file path')
     parser.add_argument('-o','--outdir',type=str,help='Output directory path')
     parser.add_argument('-t','--type', dest='phred',type=str, default='sanger', help='Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8))')
     parser.add_argument('-q','--qual',type=str, default='30', help='Minimum sequence quality threshold')
     parser.add_argument('-n','--threads', type=str, default='2', help='Number of threads allocated')
     parser.add_argument('-s','--sample',type=str, help='sample name')
     args = parser.parse_args()
     out_log = args.outdir + '/log.out'
     err_log = args.outdir + '/log.err'
     sample = [args.sample+'_1',args.sample+'_2']
     ret = run_qc(args.read1, args.read2, args.outdir,args.threads, out_log, err_log, args.phred, args.qual,sample)
