import os
import sys
import time
import datetime
import argparse
import subprocess
import config as cf
from time import gmtime, strftime

def run_varcaller(bamfile, outdir,threads, out_log, err_log, reference, vcf, intervals, db, thresh):
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('Unified Genotyper caller started at '+ start_time + '\n')
    print('Unified Genotyper started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('Unified Genotyper started at '+ start_time + '\n')
    variants = outdir+'/'+os.path.basename(bamfile)[:os.path.basename(bamfile).rfind('.')]+'.vcf'
    run_ug = subprocess.Popen(['java','-jar',cf.gatk,'-T','UnifiedGenotyper','-R',reference,'-I',bamfile,'--dbsnp',vcf,'-out_mode','EMIT_ALL_CONFIDENT_SITES','-L',intervals,'-o',variants],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    ug_status = run_ug.communicate()
    ug_exitcode = run_ug.returncode


    out_logger.write(ug_status[0] + '\n')
    err_logger.write(ug_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if ug_exitcode:
        out_logger.write('Unified Genotyper crashed at ' +end_time+' with returncode '+str(ug_exitcode)+'. Please check error log for details.\n')
        print('Unfied Genotyper crashed at ' +end_time+' with returncode '+str(ug_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Unified Genotyper crashed at '+ end_time+'\n')
        err_logger.close()
        return(ug_exitcode, variants)

    out_logger.write('Unified Genotyper completed successfully at ' +end_time)
    err_logger.write('Unified Genotyper completed successfully at ' +end_time)
    print('Unified Genotyper completed successfully at '+end_time)
    out_logger.close()
    err_logger.close()
    #End of Unified Genotyper and start of annovar

    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('Annovar started at '+ start_time + '\n')
    print('Annovar started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('Annovar started at '+ start_time + '\n')
    variants_ann = outdir+'/'+os.path.basename(bamfile)[:os.path.basename(bamfile).rfind('.')]+'_annotated.vcf'
    #add annovar command
    run_annovar = subprocess.Popen(['perl',cf.annovar+'table_annovar.pl',variants,cf.annovar+'humandb/','-buildver','hg19','-out',variants_ann,'-remove','-protocol',db,'-operation','g,r,r,f,f,f,f','-nastring','.','-vcfinput'],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    annovar_status = run_annovar.communicate()
    annovar_exitcode = run_annovar.returncode


    out_logger.write(annovar_status[0] + '\n')
    err_logger.write(annovar_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if annovar_exitcode:
        out_logger.write('Annovar crashed at ' +end_time+' with returncode '+str(annovar_exitcode)+'. Please check error log for details.\n')
        print('Annovar crashed at ' +end_time+' with returncode '+str(annovar_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('Annovar crashed at '+ end_time+'\n')
        err_logger.close()
        return(annovar_exitcode, variants_ann)

    out_logger.write('Annovar completed successfully at ' +end_time)
    err_logger.write('Annovar completed successfully at ' +end_time)
    print('Annovar completed successfully at '+end_time)
    out_logger.close()
    err_logger.close()
    return(annovar_exitcode, variants_ann)
    #End of Annovar

if __name__ == '__main__':
    parser = argparse.ArgumentParser('This module runs GATK jars to realign indels and recalibrate base quality score from bam files')
    parser.add_argument('-b','--bamfile',type=str,help='Input sam file path')
    parser.add_argument('-o','--outdir',type=str,help='Output directory path')
    parser.add_argument('-r','--reference',type=str,help='Path to reference fasta file')
    parser.add_argument('-t','--threads',type=str,help='Number of threads allocated')
    parser.add_argument('-v','--vcf',type=str,help='Path to known variants file')
    parser.add_argument('-i','--intervals',type=str, help='Path to intervals file')
    parser.add_argument('-d','--db',type=str, nargs='+', help='Databases to annotate variants with')
    parser.add_argument('-q','--thresh', type=str, nargs='+', help='Quality threshold values')
    args = parser.parse_args()
    out_log = args.outdir + '/log.out'
    err_log = args.outdir + '/log.err'
    ret = run_varcaller(args.bamfile, args.outdir, args.threads, out_log, err_log, args.reference, args.vcf, args.intervals, ','.join(args.db), args.thresh)

