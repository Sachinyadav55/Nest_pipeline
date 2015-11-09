#! /usr/bin/python3
import os
import re
import sys
import csv
import glob
import gatk
import time
import bowtie
import picard
import logging
import argparse
import pandas as pd
import trimgalore
import subprocess
import configparser
import collections
from multiprocessing import Pool
from collections import defaultdict
from SeqCohort import SeqCohort
from SeqMultiple import SeqMultiple
from SeqRapid import SeqRapid

def main():
    config = configparser.ConfigParser()
    config.read('config.cfg')
    pipeline = config['General']['pipeline']    
    if pipeline == 'dnaseq_cohort':
        dc = SeqCohort()
        output_record, file_record = dc.create_outdir()
        data = [[val,file_record[val][0],file_record[val][1],
            output_record[val]] for val in file_record]
        outdir = output_record.values()
        pool_dc = Pool(3) 
        gvcf_list = pool_dc.map(dc.preprocess, data)
        vcffile = dc.joint_genotyping(gvcf_list)
        vcffile = dc.vqsr(vcffile)

    elif pipeline == 'dnaseq':
        dc = SeqMultiple()
        output_record, file_record = dc.create_outdir()
        data = [[val,file_record[val][0],file_record[val][1],
            output_record[val]] for val in file_record]
        outdir = output_record.values()
        pool_dc = Pool(3) 
        gvcf_list = pool_dc.map(dc.preprocess, data)

    elif pipeline == 'exomeseq_cohort':
        dc = SeqCohort()
        output_record, file_record = dc.create_outdir()
        data = [[val,file_record[val][0],file_record[val][1],
            output_record[val]] for val in file_record]
        outdir = output_record.values()
        pool_dc = Pool(3) 
        gvcf_list = pool_dc.map(dc.preprocess, data)
        vcffile = dc.joint_genotyping(gvcf_list)
        vcffile = dc.vqsr(vcffile)

    elif pipeline == 'exomeseq':
        dc = SeqMultiple()
        output_record, file_record = dc.create_outdir()
        data = [[val,file_record[val][0],file_record[val][1],
            output_record[val]] for val in file_record]
        outdir = output_record.values()
        pool_dc = Pool(3) 
        gvcf_list = pool_dc.map(dc.preprocess, data)
    
    elif pipeline == 'exomeseq_rapid':
        dc = SeqRapid()
        output_record, file_record = dc.create_outdir()
        data = [[val,file_record[val][0],file_record[val][1],
            output_record[val]] for val in file_record]
        outdir = output_record.values()
        pool_dc = Pool(3) 
        fastq_list = pool_dc.map(dc.preprocess, data)
        fwdlane = defaultdict(list)
        revlane = defaultdict(list)
        for f,r in fastq_list:
            lane = os.path.basename(f).split('_L')[0]
            fwdlane[lane].append(f)
            lane = os.path.basename(r).split('_L')[0]
            revlane[lane].append(r)
        fwd = sorted([','.join(fwdlane[lanes]) for lanes in fwdlane])
        rev = sorted([','.join(revlane[lanes]) for lanes in revlane])
        base = [os.path.basename(os.path.dirname(os.path.dirname
            (os.path.dirname(val.split(',')[0])))) for val in fwd]
        outdir = [os.path.dirname(os.path.dirname(os.path.dirname(
            val.split(',')[0]))) for val in fwd]
        data = zip(base, fwd, rev, outdir)
        pool_dc = Pool(3)
        gvcf_list = pool_dc.map(dc.align, data)
    return

def configure(fwd,rev,outdir,threads,java,mem,reference,kmer,adap,window,
        lead,trail,crop,headcrop,minlen,mode,mismatch,seedlen,interval,seedex,
        reseed,ambfunc,dpad,gbar,matbonus,penrange,ambpen,gappen,refgappen,
        maplim,minins,maxins,orient,rgid,rglb,rgpl,rgpu,rgsm,rgcn,rgds,rgdt,
        rgpi,rgpg,rgpm,dup,dupscore,dupregex,duppix,maxint,minreads,mismatches,
        windows,model,lod,entropy,maxcon,maxmoves,greedy,maxreads,
        maxinmem,covariates,ics,maxcyc,mcs,bqsrpen,ddq,idq,calconf,
        emitconf,gtmode,mmq,ann,applymode,filterlevel,exome,pipeline,erc):
    
    #General config
    project = os.path.dirname(os.path.abspath(__file__))
    ref_dir = '/data/db/'
    ref_path = '{1}/{0}/bowtie/genome.fa'.format(reference,
        ref_dir)
    threads = str(int(int(threads)/3))
    config = configparser.ConfigParser()
    
    #Database paths
    data = '/data/db/{0}/gatk/'.format(reference)
    mills = '{0}Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz'.format(data)
    okg = '{0}1000G_phase1.indels.hg19.sites.vcf.gz'.format(data)
    dbsnp = '{0}dbsnp.vcf.gz'.format(data)
    hapmap = '{0}hapmap_3.3.hg19.sites.vcf.gz'.format(data)
    omni = '{0}1000G_omni2.5.hg19.sites.vcf.gz'.format(data)
    known = ','.join([mills, okg])
    knownsites = ','.join([dbsnp, mills, okg])
    knownvqsr = ','.join([hapmap, omni, okg, dbsnp])
    
    #General config
    config['General'] = {'fwd':','.join(fwd), 'rev':','.join(rev), 
        'outdir': outdir, 'threads':threads, 'java':java, 'mem':mem, 
        'reference': ref_path, 'dbsnp': dbsnp, 'pipeline':pipeline}
    
    #FastQC config
    fastqc_path = '{0}/FastQC/fastqc'.format(project)
    config['FastQC'] = {'fastqc': fastqc_path, 'kmer':kmer}

    #Trimmomatic config
    trimmomatic = '{0}/Trimmomatic-0.33/trimmomatic-0.33.jar'.format(project)
    config['Trimmomatic'] = {'trimmomatic': trimmomatic, 'illuminaclip' : adap, 
        'window': window, 'leading': lead, 'trailing': trail, 'crop': crop, 
        'headcrop' : headcrop, 'minlen' : minlen}

    #Bowtie config
    bowtiepath = '{0}/bowtie2-2.2.5/bowtie2'.format(project)
    bowtieindex = '{1}/{0}/bowtie/genome'.format(reference,
        ref_dir)
    presets = {'5,1,0,22,S,0,2.50,end-to-end' : '--very-fast',
        '10,2,0,22,S,0,2.50,end-to-end' : '--fast',
        '15,2,0,22,S,1,1.15,end-to-end' : '--sensitive',
        '20,3,0,20,S,1,0.50,end-to-end' : '--very-sensitive',
        '5,1,0,25,S,1,2.00,local' : '--very-fast-local',
        '10,2,0,22,S,1,1.75,local' : '--fast-local',
        '15,2,0,20,S,1,0.75,local' : '--sensitive-local',
        '20,3,0,20,S,1,0.50,local' : '--very-sensitive-local'}
    try:
        runmode = presets[','.join([seedex, reseed, mismatch, seedlen,
            interval, mode])] 
        config['Bowtie'] = {'bowtie': bowtiepath, 'index': bowtieindex, 
            'preset':runmode, 'ambfunc': ambfunc, 'dpad':dpad, 'gbar': gbar, 
            'matbonus': matbonus, 'penrange': penrange, 'ambpen': ambpen, 
            'gappen': gappen, 'refgappen': refgappen, 'maplim': maplim, 
            'minins': minins, 'maxins': maxins, 'orient': '--'+orient}
    except KeyError:
        config['Bowtie'] = {'bowtie': bowtiepath, 'index': bowtieindex, 
            'mode': '--'+mode, 'mismatch': mismatch, 'seedlen':seedlen,
            'interval': interval, 'seedex': seedex, 'reseed': reseed,
            'ambfunc': ambfunc, 'dpad':dpad, 'gbar': gbar, 'matbonus': matbonus,
            'penrange': penrange, 'ambpen': ambpen, 'gappen': gappen,
            'refgappen': refgappen, 'maplim': maplim, 'minins': minins,
            'maxins': maxins, 'orient': '--'+orient}
    if maplim == '1':
        del config['Bowtie']['maplim']

    #Picard config
    picard = '{0}/picard-tools-1.139/picard.jar'.format(project)
    config['Picard'] = {'picard': picard, 'rgid':rgid, 'rglb': rglb,
        'rgpl': rgpl, 'rgpu' : rgpu, 'rgsm': rgsm, 'rgcn': rgcn, 'rgds': rgds,
        'rgdt':rgdt, 'rgpi' : rgpi, 'rgpg' : rgpg, 'rgpm': rgpm, 'dup': dup,
        'dupscore': dupscore, 'dupregex': dupregex, 'duppix': duppix}

    #GATK config
    gatk = '{0}/GenomeAnalysisTK/GenomeAnalysisTK.jar'.format(project)
    if mode == 'end-to-end':
        mmq = 0
    config['GATK'] = {'gatk': gatk, 'maxins': maxint, 'minreads': minreads,
        'mismatch': mismatches, 'window': windows, 'known' : known,
        'model': model, 'lod': lod, 'entropy': entropy, 'maxcon': maxcon,
        'maxmoves': maxmoves, 'greedy': greedy, 'maxreads': maxreads,
        'mmq': mmq, 'maxinmem': maxinmem, 'covariates': covariates,
        'knownsites': knownsites, 'ics': ics, 'maxcyc': maxcyc, 'mcs': mcs,
        'bqsrpen': bqsrpen, 'ddq': ddq, 'idq': idq, 'calconf': calconf,
        'emitconf': emitconf, 'gtmode': gtmode, 'knownvqsr': knownvqsr, 
        'ann': ann, 'applymode': applymode, 'filterlevel': filterlevel,
        'exome':exome, 'erc':erc}

    #Write config file
    with open('config.cfg','w') as configfile:
        config.write(configfile)
    return(os.path.abspath(project+'/config.cfg'))

if __name__ == '__main__' :
    nester = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    general = nester.add_argument_group('General')
    general.add_argument('-i', '--indir', type=str, help='Input directory path.')
    general.add_argument('-c', '--config', type=str, 
        default=os.path.dirname(os.path.abspath(__file__))+'.config.cfg', 
        help='Config file path.')
    general.add_argument('--outdir', type=str, default=None,
        help='Output directory path.')
    general.add_argument('--thread', type=str, default='6', help='Number of threads.')
    general.add_argument('--java', type=str, default='/usr/bin/java', help='Java path.')
    general.add_argument('--mem', type=str, default='4g', help='java memory limit.')
    general.add_argument('--reference', type=str, default='hg19',
        help='Genome build to use for analysis.')
    capture = '/data/db/hg19/exome/nexterarapidcapture_exome_targetedregions_v1.2.bed'
    general.add_argument('--pipeline', type=str, default='exomeseq_cohort',
        choices=['exomeseq_cohort', 'dnaseq_cohort', 'dnaseq', 'exomeseq', 
        'exomeseq_rapid'], help='Pipeline to run.')
    general.add_argument('--exome', type=str, default=capture,
        help='Exome capture kit interval file.')
    fastqc = nester.add_argument_group('FastQC')
    fastqc.add_argument('--fastqckmer', type=str, default='7', 
        help='FastQC kmer lenght to analyse.')
    trim = nester.add_argument_group('Trimmomatic')
    trim.add_argument('--window', type=str, default='4:15',
        help='Trimmomatic sliding window and quality threshold.')
    trim.add_argument('--cropl', type=str, default='0', 
        help='Length of bases at 5\' end to be dicarded if below quality threshold.')
    trim.add_argument('--cropt', type=str, default='0', 
        help='Length of bases at 3\' end to be discarded if below quality therhold.')
    trim.add_argument('--headcrop', type=str, default='0',
        help='Clip n bases from the 5\' end')
    trim.add_argument('--minlen', type=str, default='36',
        help='Minimum read length to retain read')
    trim.add_argument('--adapters', type=str, default='adapters.txt',
        help='File containing adapter seqeunces.')
    bowties = nester.add_argument_group('Bowtie')
    bowties.add_argument('--mode', type=str, default='local',
        help='Bowtie2 alignment mode.', choices=['local', 'global'])
    bowties.add_argument('--mismatch', type=str, default='0',
        help='Mismatches allowed')
    bowties.add_argument('--seedlen', type=str, default='20', help='Seed length')
    bowties.add_argument('--interval', type=str, default='S,1,0.50',
        help='Seed extension function')
    bowties.add_argument('--seedx', type=str, default='20', help='Seed extension')
    bowties.add_argument('--reseed', type=str, default='3', help='Reseed value')
    bowties.add_argument('--ambfunc', type=str, default='L,0,0.15', 
        help='A function to govern the maximum number of ambiguous characters.')
    bowties.add_argument('--dpad', type=str, default='15', 
        help='"Pads" dynamic programming problems, by n colmuns to allow gaps.')
    bowties.add_argument('--gbar', type=str, default='4',
        help='Disallow gaps n bases from start or end of read.')
    bowties.add_argument('--matbonus', type=str, default='2',
        help='Set match bonus.')
    bowties.add_argument('--penrange', type=str, default='6,2',
        help='Set maximum and minimum mismatch penalties.')
    bowties.add_argument('--ambpen', type=str, default='1',
        help='Set penality for positions in read and reference with N\'s')
    bowties.add_argument('--gappen', type=str, default='5,3',
        help='Set read gap open and extend penalties.')
    bowties.add_argument('--refgappen', type=str, default='5,3',
        help='Set reference gap open and extend penalties.')
    bowties.add_argument('--maplim', type=str, default='1',
        help='Number of distinct alignments allowed for a read.')
    bowties.add_argument('--minins', type=str, default='0', 
        help='Minimum fragment length for valid paired end alignment.')
    bowties.add_argument('--maxins', type=str, default='500',
        help='Maximum fragment length for valid paired end alignment.')
    bowties.add_argument('--orient', type=str, default='fr',
        help='Paired end read orientation.')
    arg = nester.add_argument_group('Add Or Replace Read Group')
    arg.add_argument('--rgid', type=str, default='null',
        help='Read group ID.')
    arg.add_argument('--rglb', type=str, default='Exome_Seq',
        help='Read group library.')
    arg.add_argument('--rgpl', type=str, default='Illumina',
        help='Read group platform.')
    arg.add_argument('--rgpu', type=str, default='HiSeq2500',
        help='Read group platform unit.')
    arg.add_argument('--rgsm', type=str, default='null',
        help='Read group sample name.')
    arg.add_argument('--rgcn', type=str, default='VannbergLab',
        help='Read group sequencing center name.')
    arg.add_argument('--rgds', type=str, default='null',
        help='Read group description.')
    arg.add_argument('--rgdt', type=str, default='null',
        help='Read group run date.')
    arg.add_argument('--rgpi', type=str, default='null',
        help='Read group predicted insert size.')
    arg.add_argument('--rgpg', type=str, default='null',
        help='Read group program group.')
    arg.add_argument('--rgpm', type=str, default='null',
        help='Read group platform model')
    mark = nester.add_argument_group('Mark Duplicates')
    mark.add_argument('--dup', type=str, default='false',
        help='If true will not write duplicates in to bam file.')
    mark.add_argument('--dupscore', type=str, default='SUM_OF_BASE_QUALITIES',
        help='The scoring strategy for choosing non-duplicate reads.')
    mark.add_argument('--dupregex', type=str,
        default='[a-zA-Z0-9]+:[0-9]:([0-9]+):(0-9]+).*', 
        help='Regex expression that can be used to parse read names.')
    mark.add_argument('--duppix', type=str, default='100',
        help='Maximum offset between two duplicate clusters to consider\
        them optical duplicates.')
    indel = nester.add_argument_group('Indel realignment')
    indel.add_argument('--maxint', type=str, default='500', 
        help='Maximum interval size.')
    indel.add_argument('--minreads', type=str, default='4',
        help='Minimum reads at a locus to enable entropy calculation.')
    indel.add_argument('--mismatches', type=str, default='0.0',
        help='Fraction of base qualities needing to mismatch for a position\
        to have high entropy.')
    indel.add_argument('--windows', type=str, default='10',
        help='Window size for calculating entropy or SNP clusters.')
    indel.add_argument('--model', type=str, default='USE_READS',
        help='Determines how to compute the possible alternate consenses.')
    indel.add_argument('--lod', type=str, default='5.0',
        help='LOD threshold above which the cleaner will clean.')
    indel.add_argument('--entropy', type=str, default='0.15',
        help='Percentage of a mismatches at a locus to be considered having\
        high entropy.')
    indel.add_argument('--maxcon', type=str, default='30',
        help='Max alternate consensuses to try.')
    indel.add_argument('--maxmoves', type=str, default='200',
        help='Maximum positional move in basepairs that a read can be adjusted\
        during realignment.')
    indel.add_argument('--greedy', type=str, default='120',
        help='Max reads used for finding the alternate consensuses.')
    indel.add_argument('--maxreads', type=str, default='20000',
        help='Max reads allowed at an interval for realignment.')
    indel.add_argument('--maxinmem', type=str, default='150000', 
        help='Max reads allowed to be kept in memory at a time by the SAMFileWriter.')
    base = nester.add_argument_group('Base Quality Score Recalibration')
    base.add_argument('--covariantes', type=str, nargs='+', 
        default='ReadGroupCovariate,QualityScoreCovariate,CycleCovariate,ContextCovariate',
        help='Covariatesi to be used in the recalibration.')
    base.add_argument('--ics', type=str, default='3',
        help='Size of k-mer context to be used for base indels.')
    base.add_argument('--maxcyc', type=str, default='500', 
        help='The maximum cycle value permitted for cycle covariates.')
    base.add_argument('--mcs', type=str, default='2',
        help='Size of the k-mer context to be used for base mismatches.')
    base.add_argument('--bqsrgpen', type=str, default='40.0',
        help='BQSR BAQ gap open penalty.')
    base.add_argument('--ddq', type=str, default='45',
        help='Default quality for the base deletions covariate.')
    base.add_argument('--idq', type=str, default='45',
        help='Default quality for base insertion covariate.')
    haplo = nester.add_argument_group('Haplotype Caller')
    haplo.add_argument('--calconf', type=str, default='30.0', 
        help='The minimum phred scaled confidence threshold for calling variants.')
    haplo.add_argument('--emitconf', type=str, default='10.0',
        help='The minimum phred scaled confidence threshold for emitting variants.')
    haplo.add_argument('--gtmode', type=str, default='DISCOVERY',
        help='Specify how to determine the alternate allele.')
    haplo.add_argument('--mmq', type=str, default='20', 
        help='Minimum read mapping quality for variant calling.')
    haplo.add_argument('--erc', type=str, default='NONE',
        help='Mode for emitting reference confiedence scores.')
    vqsr = nester.add_argument_group('Variant Quality Score Recalibration')
    vqsr.add_argument('--ann', type=str, 
        default='QD,MQ,MQRankSum,ReadPosRankSum,FS,SOR',
        help='Annotations used for vairant score recalibration')
    vqsr.add_argument('--applymode', type=str, default='BOTH',
        help='VQSR apply mode')
    vqsr.add_argument('--filterlevel', type=str, default='99.0',
        help='Filter level')
    args = nester.parse_args()
    if args.outdir == None:
        args.outdir = args.indir + 'Outputs'
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    if args.pipeline == 'Exome_Seq' and args.exome == 'null':
        logger.error('Exome capture intervals not provided')
        sys.exit()
    read = re.compile('_[Rr]*1')
    basename = re.compile('_[Rr]*[12]')
    fastq = re.compile('fastq|fq')
    input_files = glob.glob(os.path.abspath(args.indir+'*.fastq*'))
    if len(input_files) == 0:
        input_files = glob.glob(os.path.abspath(args.indir+'/*/*.fastq*'))
    input_dict = defaultdict(list)
    for lines in input_files:
        base = basename.split(fastq.split(os.path.basename(lines))[0])[0]
        input_dict[base].append(lines)
    fwd = list()
    rev = list()
    for lines in input_dict:
        if read.search(input_dict[lines][0]):
            fwd.append(input_dict[lines][0])
            rev.append(input_dict[lines][1])
        else:
            fwd.append(input_dict[lines][1])
            rev.append(input_dict[lines][0])
            
    if args.rgdt == 'null':
            args.rgdt = time.strftime('%Y-%m-%d')
    configure(fwd,rev,args.outdir,args.thread,args.java,args.mem,args.reference,
            args.fastqckmer, args.adapters, args.window, args.cropl, args.cropt, 
            args.headcrop, args.headcrop, args.minlen, args.mode, args.mismatch,
            args.seedlen, args.interval, args.seedx, args.reseed, args.ambfunc,
            args.dpad, args.gbar, args.matbonus, args.penrange, args.ambpen,
            args.gappen, args.refgappen, args.maplim, args.minins, args.maxins, 
            args.orient, args.rgid, args.rglb, args.rgpl, args.rgpu, args.rgid, 
            args.rgcn, args.rglb, args.rgdt, args.rgpi, args.rgpg, args.rgpm, 
            args.dup, args.dupscore, args.dupregex, args.duppix,args.maxint, 
            args.minreads, args.mismatches, args.windows, 
            args.model, args.lod, args.entropy, args.maxcon, args.maxmoves,
            args.greedy, args.maxreads,args.maxinmem, args.covariantes,
            args.ics, args.maxcyc, args.mcs, args.bqsrgpen,args.ddq,
            args.idq, args.calconf, args.emitconf, args.gtmode, args.mmq, 
            args.ann, args.applymode, args.filterlevel, args.exome,
            args.pipeline, args.erc)
    main()
