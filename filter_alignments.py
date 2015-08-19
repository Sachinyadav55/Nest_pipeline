import os
import sys
import time
import datetime
import argparse
import subprocess
import config as cf
from time import gmtime, strftime

def run_picard(samfile, outdir, out_log, err_log, platform, sample, lane):
    #Start of sorting of bam files
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    out_logger.write('SortSam started at '+ start_time + '\n')
    print('SortSam started at '+ start_time )
    err_logger = open(err_log, 'a')
    err_logger.write('SortSam started at '+ start_time + '\n')
    bamfile = outdir + '/' + os.path.basename(samfile)[:os.path.basename(samfile).rfind('.')] + '_SO.bam'
    run_sortsam = subprocess.Popen(['java','-jar',cf.picard+'SortSam.jar','I='+samfile,'O='+bamfile,'SO=coordinate','CREATE_INDEX=true'],stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    sortsam_status = run_sortsam.communicate()
    sortsam_exitcode = run_sortsam.returncode

    out_logger.write(sortsam_status[0] + '\n')
    err_logger.write(sortsam_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if sortsam_exitcode:
        out_logger.write('SortSam crashed at ' +end_time+' with returncode '+str(sortsam_exitcode)+'. Please check error log for details.\n')
        print('SortSam crashed at ' +end_time+' with returncode '+str(sortsam_exitcode)+'. Please check error log for details.')
        out_logger.close()
        err_logger.write('SortSam crashed at '+ end_time+'\n')
        err_logger.close()
        return(sortsam_exitcode, bamfile)
    out_logger.write('SortSam completed successfully at ' +end_time+'\n')
    err_logger.write('SortSam completed successfully at ' +end_time+'\n')
    print('SortSam completed successfully at '+end_time)
    out_logger.close()
    err_logger.close()

    #End of sorting and start of read group addition
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    err_logger = open(err_log, 'a')
    out_logger.write('AddOrReplaceReadGroup started at '+start_time+'\n')
    print('AddOrReplaceReadGroup started at '+start_time)
    err_logger.write('AddOrReplaceReadGroup started at '+start_time+'\n')
    
    run_adrg = subprocess.Popen(['java','-jar',cf.picard +'AddOrReplaceReadGroups.jar','I='+bamfile,'O='+outdir+'/'+os.path.basename(samfile)[:os.path.basename(samfile).rfind('.')] + '_RG.bam','VALIDATION_STRINGENCY=SILENT','RGID='+sample[:-1],'RGLB='+lane,'RGPL='+platform,'RGPU=bc','RGSM='+sample[:-1],'CREATE_INDEX=true'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    adrg_status = run_adrg.communicate()
    adrg_exitcode = run_adrg.returncode

    out_logger.write(adrg_status[0] + '\n')
    err_logger.write(adrg_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if adrg_exitcode:
        out_logger.write('AddOrReplaceReadGroup crashed at '+end_time+' with returncode '+str(adrg_exitcode)+'.Please check error log for details. \n')
        print('AddOrReplaceReadGroup crashed at '+end_time+' with returncode '+ str(adrg_exitcode)+' .Please check error log for details.\n')
        err_logger.write('AddOrReplaceReadGroup crashed at '+end_time+'\n')
        out_logger.close()
        err_logger.close()
        return(adrg_exitcode,outdir+'/'+os.path.basename(samfile)[:os.path.basename(samfile).rfind('.')] + '_RG.bam')
    out_logger.write('AddOrReplaceReadGroup completed successfully at ' +end_time+'\n')
    err_logger.write('AddOrReplaceReadGroup completed successfully at ' +end_time+'\n')
    print('AddOrReplaceReadGroup completed successfully at '+end_time+'\n')
    out_logger.close()
    err_logger.close()

 
    #End of read group addition and start duplicate removal
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    out_logger = open(out_log, 'a')
    err_logger = open(err_log, 'a')
    out_logger.write('Mark duplicates started at '+start_time+'\n')
    print('Mark duplicates started at '+start_time)
    err_logger.write('Mark duplicates started at '+start_time+'\n')

    run_rmdup = subprocess.Popen(['java','-jar',cf.picard+'MarkDuplicates.jar','I='+outdir+'/'+os.path.basename(samfile)[:os.path.basename(samfile).rfind('.')] + '_RG.bam','O='+outdir+'/'+os.path.basename(samfile)[:os.path.basename(samfile).rfind('.')] + '_RD.bam','M='+os.path.basename(samfile)[:os.path.basename(samfile).rfind('.')]+'.metrics','REMOVE_DUPLICATES=true','AS=true','CREATE_INDEX=true'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    rmdup_status = run_rmdup.communicate()
    rmdup_exitcode = run_rmdup.returncode

    out_logger.write(rmdup_status[0] + '\n')
    err_logger.write(rmdup_status[1] + '\n')
    end_time = strftime('%Y-%m-%d %H:%M:%S')
    if adrg_exitcode:
        out_logger.write('MarkDuplicates crashed at '+end_time+' with returncode '+str(adrg_exitcode)+'.Please check error log for details. \n')
        print('MarkDuplicates crashed at '+end_time+' with returncode '+ str(adrg_exitcode)+' .Please check error log for details.\n')
        err_logger.write('MarkDuplicates crashed at '+end_time+'\n')
        out_logger.close()
        err_logger.close()
        return(rmdup_exitcode, outdir+'/'+os.path.basename(samfile)[:os.path.basename(samfile).rfind('.')] + '_RD.bam')
    out_logger.write('Mark Duplicates completed successfully at ' +end_time+'\n')
    err_logger.write('Mark Duplicates completed successfully at ' +end_time+'\n')
    print('Mark Duplicates completed successfully at '+end_time)
    out_logger.close()
    err_logger.close()

    return(rmdup_exitcode, outdir+'/'+os.path.basename(samfile)[:os.path.basename(samfile).rfind('.')] + '_RD.bam')

if __name__ == '__main__':
    parser = argparse.ArgumentParser('This module runs picards jars to sort and remove duplicate records from bam files')
    parser.add_argument('-s','--sam',type=str,help='Input sam file path')
    parser.add_argument('-o','--outdir',type=str,help='Output directory path')
    parser.add_argument('-p','--platform',type=str,help='Platform of sequencing')
    parser.add_argument('-n','--sample',type=str,help='Sample name')
    parser.add_argument('-l','--lane',type=str,help='Lane information',default='l001')
    args = parser.parse_args()
    out_log = args.outdir + '/log.out'
    err_log = args.outdir + '/log.err'
    run_picard(args.sam, args.outdir, out_log, err_log, args.platform, args.sample, args.lane)
