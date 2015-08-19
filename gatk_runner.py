import os
import sys
import time
import datetime
import argparse
import subprocess
import config as cf
from time import gmtime, strftime

def run_gatk(bamfile, outdir,threads, out_log, err_log, reference, vcf):
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('Realigner target creator started at '+ start_time + '\n')
    print('Realigner target creator started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('Realigner target creator started at '+ start_time + '\n')
    intervals = outdir+'/'+os.path.basename(bamfile)[:os.path.basename(bamfile).rfind('.')]+'.intervals'
    run_rtc = subprocess.Popen(['java','-jar',cf.gatk,'-T','RealignerTargetCreator','-nt',threads,'-R',reference,'-I',bamfile, '-o', intervals],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    rtc_status = run_rtc.communicate()
    rtc_exitcode = run_rtc.returncode


    out_logger.write(rtc_status[0] + '\n')
    err_logger.write(rtc_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if rtc_exitcode:
        out_logger.write('Realigner target creator crashed at ' +end_time+' with returncode '+str(rtc_exitcode)+'. Please check error log for details.\n')
        print('Realigner target creator crashed at ' +end_time+' with returncode '+str(rtc_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Indel aligner crashed at '+ end_time+'\n')
        err_logger.close()
        return(rtc_exitcode, None,intervals)

    out_logger.write('Realigner target creator completed successfully at ' +end_time+'\n')
    err_logger.write('Realigner target creator completed successfully at ' +end_time+'\n')
    print('Realigner target creator completed successfully at '+end_time+'\n')
    out_logger.close()
    err_logger.close()
    #End of realigner target creator

    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('Indel realigner started at '+ start_time + '\n')
    print('Indel realigner started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('Indel realigner started at '+ start_time + '\n')
    run_ir = subprocess.Popen(['java','-jar',cf.gatk,'-T','IndelRealigner','-R',reference,'-I',bamfile,'-targetIntervals',intervals,'-o', bamfile[:bamfile.rfind('.')]+'_ir.bam'],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    ir_status = run_ir.communicate()
    ir_exitcode = run_ir.returncode


    out_logger.write(ir_status[0] + '\n')
    err_logger.write(ir_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if ir_exitcode:
        out_logger.write('Indel realigner crashed at ' +end_time+' with returncode '+str(ir_exitcode)+'. Please check error log for details.\n')
        print('Indel realigner crashed at ' +end_time+' with returncode '+str(ir_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Indel realigner crashed at '+ end_time+'\n')
        err_logger.close()
        return(ir_exitcode, bamfile[:bamfile.rfind('.')]+'_ir.bam',None)

    out_logger.write('Indel realigner completed successfully at ' +end_time+'\n')
    err_logger.write('Indel realigner completed completed successfully at ' +end_time+'\n')
    print('Indel realigner completed successfully at '+end_time+'\n')
    out_logger.close()
    err_logger.close()
    bamfile = bamfile[:bamfile.rfind('.')]+'_ir.bam'
    #End of indel realigner
    #Start of base recalibrator
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('Base quality score recalibrator started at '+ start_time + '\n')
    print('Base quality score recalibrator started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('Base quality score recalibrator started at '+ start_time + '\n')
    run_br = subprocess.Popen(['java','-jar',cf.gatk,'-T','BaseRecalibrator','-R',reference,'-I',bamfile,'--knownSites',vcf,'-o', bamfile[:bamfile.rfind('.')]+'.table'],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    br_status = run_br.communicate()
    br_exitcode = run_br.returncode


    out_logger.write(br_status[0] + '\n')
    err_logger.write(br_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if br_exitcode:
        out_logger.write('Base recalibrator crashed at ' +end_time+' with returncode '+str(br_exitcode)+'. Please check error log for details.\n')
        print('Base recalibrator crashed at ' +end_time+' with returncode '+str(br_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Base recalibrator crashed at '+ end_time+'\n')
        err_logger.close()
        return(br_exitcode, bamfile,None)

    out_logger.write('Base recalibrator completed successfully at ' +end_time+'\n')
    err_logger.write('Base recalibrator completed successfully at ' +end_time+'\n')
    print('Base recalibrator completed successfully at '+end_time+'\n')
    out_logger.close()
    err_logger.close()
    #End of base recalibrator 
    #Start of print reads
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('Print reads started at '+ start_time + '\n')
    print('Print reads started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('Print reads started at '+ start_time + '\n')
    run_pr = subprocess.Popen(['java','-jar',cf.gatk,'-T','PrintReads','-R',reference,'-I',bamfile,'--BQSR',bamfile[:bamfile.rfind('.')]+'.table','-o', bamfile[:bamfile.rfind('_')]+'_recal.bam'],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    pr_status = run_pr.communicate()
    pr_exitcode = run_pr.returncode


    out_logger.write(pr_status[0] + '\n')
    err_logger.write(pr_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if pr_exitcode:
        out_logger.write('Base recalibrator crashed at ' +end_time+' with returncode '+str(pr_exitcode)+'. Please check error log for details.\n')
        print('Base recalibrator crashed at ' +end_time+' with returncode '+str(pr_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Base recalibrator crashed at '+ end_time+'\n')
        err_logger.close()
        return(pr_exitcode, bamfile[:bamfile.rfind('_')]+'_recal.bam',None)

    out_logger.write('Print Reads completed successfully at ' +end_time+'\n')
    err_logger.write('Print Reads completed successfully at ' +end_time+'\n')
    print('Print Reads completed successfully at '+end_time+'\n')
    out_logger.close()
    err_logger.close()

    return(pr_exitcode,bamfile[:bamfile.rfind('_')]+'_recal.bam',intervals)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('This module runs GATK jars to realign indels and recalibrate base quality score from bam files')
    parser.add_argument('-b','--bamfile',type=str,help='Input sam file path')
    parser.add_argument('-o','--outdir',type=str,help='Output directory path')
    parser.add_argument('-r','--reference',type=str,help='Path to reference fasta file')
    parser.add_argument('-t','--threads',type=str,help='Number of threads allocated')
    parser.add_argument('-v','--vcf',type=str,help='Path to known variants file')
    args = parser.parse_args()
    out_log = args.outdir + '/log.out'
    err_log = args.outdir + '/log.err'
    ret = run_gatk(args.bamfile, args.outdir, args.threads, out_log, err_log, args.reference, args.vcf)

