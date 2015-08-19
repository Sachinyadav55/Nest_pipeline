import os
import sys
import time
import datetime
import argparse
import subprocess
import filter_fastq
import align_reads
import filter_alignments
import gatk_runner
import multiprocessing
import variant_detection
import cosmos
import config as cf
from time import gmtime, strftime
from multiprocessing import Pool
from itertools import repeat

def run_seq(arguments):
    fastq1 = arguments[0]
    fastq2 = arguments[1]
    outdir = arguments[2]
    phred = arguments[3]
    platform = arguments[4]
    fqual = arguments[5]
    threads = arguments[6]
    genome = arguments[7]
    sample = arguments[8]
    var_qual = arguments[9]
    var_maf = arguments[10]
    lane = arguments[11]
    alignment = arguments[12]
    out_log = outdir + '/out.log'
    err_log = outdir + '/err.log'
    sample_list = [sample+'_1',sample+'_2']
    align_mode = {'local':'--very-sensitive-local','global':'--very-sensitive'}
    start_time = strftime('%Y-%m-%d %H:%M:%S')
    print('Seq analyze pipeline started at '+ start_time)
    qual_convert = {'phred64':'illumina', 'phred33':'sanger'}
    reference = {'hg19':[cf.hg19+'/Sequence/WholeGenomeFasta/genome.fa',cf.dbsnp+'/dbsnp_138.hg19.vcf','refGene,cytoBand,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp138,ljb23_all',cf.hg19+'/Sequence/Bowtie2Index/genome']}
    print(reference[genome])
    out_logger = open(out_log,'a')
    err_logger = open(err_log,'a')
    out_logger.write('Seq analyze pipeline started at '+ start_time +'\n')
    err_logger.write('Seq analyze pipeline started at '+ start_time +'\n')
    out_logger.close()
    err_logger.close()
    #Run fastq qc and fastq trimming
    fastq_filterer = filter_fastq.run_qc(fastq1, fastq2, outdir,threads, out_log, err_log, qual_convert[phred], fqual, sample_list)
    if fastq_filterer[0]:
        end_time = strftime('%Y-%m-%s %H:%M:%S')
	out_logger = open(out_log, 'a')
        err_logger = open(err_log, 'a')
        print('Seq analyze crashed at '+ end_time +'. Please check error log for details.')
        out_logger.write('Seq analyze crashed at '+ end_time +'. Please check error log for details.\n')
        err_logger.write('Seq analyze crashed at '+ end_time +'.\n')
        out_logger.close()
        err_logger.close()
	sys.exit()
    
    else:
        fastq1 = fastq_filterer[1]
        fastq2 = fastq_filterer[2]
    #Run alignment against reference genome
    fastq_aligner = align_reads.read_aligner(fastq1,fastq2,reference[genome][3],threads,outdir,out_log,err_log, phred,sample,align_mode[alignment])
    if fastq_aligner[0]:
        end_time = strftime('%Y-%m-%s %H:%M:%S')
        out_logger = open(out_log, 'a')
        err_logger = open(err_log, 'a')
        print('Seq analyze crashed at '+ end_time +'. Please check error log for details.')
        out_logger.write('Seq analyze crashed at '+ end_time +'. Please check error log for details.\n')
        err_logger.write('Seq analyze crashed at '+ end_time +'.\n')
        out_logger.close()
        err_logger.close()
        sys.exit()

    else:
	samfile = fastq_aligner[1]
    #Run picard tools to sort and add read group information to alignment files
    bam_processer = filter_alignments.run_picard(samfile, outdir, out_log, err_log, platform, sample,lane)
    if bam_processer[0]:
        end_time = strftime('%Y-%m-%s %H:%M:%S')
        out_logger = open(out_log, 'a')
        err_logger = open(err_log, 'a')
        print('Seq analyze crashed at '+ end_time +'. Please check error log for details.')
        out_logger.write('Seq analyze crashed at '+ end_time +'. Please check error log for details.\n')
        err_logger.write('Seq analyze crashed at '+ end_time +'.\n')
        out_logger.close()
        err_logger.close()
        sys.exit()
    else:
        bamfile = bam_processer[1]

    #Run gatk for inde realignment and base quality score recalibration
    gatk_pipeline =  gatk_runner.run_gatk(bamfile, outdir,threads, out_log, err_log, reference[genome][0], reference[genome][1])
    if gatk_pipeline[0]:
        end_time = strftime('%Y-%m-%s %H:%M:%S')
        out_logger = open(out_log, 'a')
        err_logger = open(err_log, 'a')
        print('Seq analyze crashed at '+ end_time +'. Please check error log for details.')
        out_logger.write('Seq analyze crashed at '+ end_time +'. Please check error log for details.\n')
        err_logger.write('Seq analyze crashed at '+ end_time +'.\n')
        out_logger.close()
        err_logger.close()
        sys.exit()
    else:
        bamfile = gatk_pipeline[1]
        intervals = gatk_pipeline[2]

    #Run variant caller and variant annotater
    var_param = '--maf '+str(var_maf) + ' --qual '+str(var_qual)
    var_caller = variant_detection.run_varcaller(bamfile, outdir,threads, out_log, err_log, reference[genome][0], reference[genome][1], intervals, reference[genome][2], var_param)
    if var_caller[0]:
        end_time = strftime('%Y-%m-%s %H:%M:%S')
        out_logger = open(out_log, 'a')
        err_logger = open(err_log, 'a')
        print('Seq analyze crashed at '+ end_time +'. Please check error log for details.')
        out_logger.write('Seq analyze crashed at '+ end_time +'. Please check error log for details.\n')
        err_logger.write('Seq analyze crashed at '+ end_time +'.\n')
        out_logger.close()
        err_logger.close()
        sys.exit()
    else:
        end_time = strftime('%Y-%m-%s %H:%M:%S')
        out_logger = open(out_log, 'a')
        err_logger = open(err_log, 'a')
        print('Seq analyze completed at '+ end_time +'. Annotated variant file : '+ var_caller[1])
        out_logger.write('Seq analyze completed successfully at '+ end_time +'.\n')
        err_logger.write('Seq analyze completed successfully at '+ end_time +'.\n')
        out_logger.close()
        err_logger.close()

    cosmic_csv = csv.writer(open(var_caller[:var_caller.rfind('.')]+'_onco_mutations.csv','w'),delimiter='\t')
    cosmic_csv.writerow(['CHROM','POSITION','REF','ALT','ONCOGENE'])
    for chrom in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']:
        onco_genes = cosmos.cosmic_reader(chrom)
        vcf_records = cosmos.vcf_reader(var_caller,chrom)
        cosmic_calls = cosmos.merge(vcf_records, onco_genes)
        for lines in cosmic_calls:
            cosmic_csv.writerow(lines)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Seq analyze takes in fastq file and run parameters and gives annotated vcf files')
    parser.add_argument('-f','--forward',nargs='+',type=str,help='Path to forward read')
    parser.add_argument('-r','--reverse',nargs='+',type=str,help='Path to reverse read')
    parser.add_argument('-o','--outdir',type=str, help='Path to output directory')
    parser.add_argument('-p','--phred', type=str, help='Type of quality values (phred33 or phred64)')
    parser.add_argument('-l','--platform', type=str, help='Type of platform used for sequencing (illumina,iontorrent,pacbio etc)')
    parser.add_argument('--lane',nargs='+',type=str, help='Lane information', default='l001')
    parser.add_argument('-q','--qual', type=str, help='Quality threshold for fastq trimming', default='30')
    parser.add_argument('-t','--thread', type=str, help='Number of threads')
    parser.add_argument('-g','--genome',type=str, help='Build of genome to be used for alignment')
    parser.add_argument('-s','--sample',nargs='+',type=str, help='Sample name')
    parser.add_argument('-v','--var_qual',type=str, help='Quality threshold for variant called')
    parser.add_argument('-m','--maf',type=str, help='Threshold for MAF of variant detected')
    parser.add_argument('-a','--alignment_mode',type=str, help='Alignment mode',choices=['global','local'],default='local') 
    args = parser.parse_args()
    print(args)
    outdir = [args.outdir + '/' + sam+lane +'/' for sam,lane in zip(args.sample,args.lane)]
    [os.mkdir(output_dir) for output_dir in outdir if not os.path.exists(os.path.abspath(output_dir))]
    cpu_count = multiprocessing.cpu_count() -1
    if len(args.sample) > 1:
        threads = cpu_count / int(args.thread)
    else:
        threads = 1
    pool = Pool(processes=threads)
    pool.map(run_seq,zip(args.forward, args.reverse, outdir, repeat(args.phred), repeat(args.platform), repeat(args.qual), repeat(args.thread), repeat(args.genome), args.sample, repeat(args.var_qual), repeat(args.maf), args.lane, repeat(args.alignment_mode)))
