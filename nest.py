import os
import re
import sys
import csv
import glob
import gatk
import time
import bowtie
import config
import picard
import logging
import argparse
import trimgalore
import subprocess
import configparser
from multiprocessing import Pool
from collections import defaultdict

class DnaseqCohort:

    read = re.compile('_[Rr]*1')
    output_record = defaultdict()
    file_record = defaultdict()

    def __init__(self):
        config = configparser.ConfigParser()
        config.read('config.cfg')

        self.read1 = config['General']['fwd'].split(',')
        self.read2 = config['General']['rev'].split(',')
        if len(self.read2) != 0:
            self.samples = [read.split(os.path.basename(val))[0] \
                            for val in self.read1]
        elif config.samples == None and len(self.read2) == 0:
            self.samples = [os.path.splitext(os.path.basename(val))[0] \
                            for val in self.read1]
        self.outdir = os.path.abspath(config['General']['outdir'])
        self.reference = config['General']['reference']
        
    def create_outdir(self):

        self.outpaths = [self.outdir + '/' + val for val in \
                        self.samples]
        output_record = dict()
        file_record = dict()
        for paths,fwd,rev,name in zip(self.outpaths,self.read1, self.read2, self.samples):
            if not os.path.exists(paths):
                os.mkdir(paths)
            output_record[name] = paths
            file_record[name] = [fwd,rev]
        return (output_record, file_record)
    
    def preprocess(self, preargs):
        self.base = preargs[0]
        self.fwd = preargs[1]
        self.rev = preargs[2]
        self.datadir = preargs[3]

        #Run quality control ananlysis and quality based trimming
        qcdir = os.path.abspath(self.datadir) + '/qc'
        if not os.path.exists(qcdir):
            os.mkdir(qcdir)
        self.fwd, self.rev = trimgalore.run_qc(self.fwd, self.rev, qcdir, self.base)

        #Run alignment
        self.bamfile = bowtie.bowtie(self.fwd, self.rev, self.datadir, self.base)
        
        #Run picard tools 
        self.bamfile = picard.addreadgroup(self.bamfile, self.base)
        self.bamfile, self.metrics = picard.markdup(self.bamfile)

        #Run GATK realignment and quality recalibration
        self.interval = gatk.target(self.bamfile)
        self.bamfile = gatk.realigner(self.bamfile, self.interval)
        self.table = gatk.baserecal(self.bamfile)
        print(self.bamfile, self.table)
        self.bamfile = gatk.printreads(self.bamfile, self.table)
        self.gvcf = gatk.haplocaller(self.bamfile)
        return(self.gvcf)

    def joint_genotyping(self, gvcf_list):
        self.gvcf_list = gvcf_list
        self.vcffile = gatk.joint_genotyper(gvcf_list, self.outdir)
        return(self.vcffile)

def main():
    config = configparser.ConfigParser()
    config.read('config.cfg')
    pipeline = config['General']['pipeline']    
    if pipeline == 'dnaseq_cohort':
        dc = DnaseqCohort()
        output_record, file_record = dc.create_outdir()
        data = [[val,file_record[val][0],file_record[val][1],output_record[val]] for val in file_record]
        #rev = [file_record[val][1] for val in file_record]
        outdir = output_record.values()
        pool_dc = Pool(2) #int(config['General']['threads']))
        gvcf_list = pool_dc.map(dc.preprocess, data)
        vcffile = dc.joint_genotyping(gvcf_list)
    return

def configure(fwd,rev,outdir,threads,java,mem,reference,kmer,adap,window,
        lead,trail,crop,headcrop,minlen,mode,mismatch,seedlen,interval,seedex,
        reseed,ambfunc,dpad,gbar,matbonus,penrange,ambpen,gappen,refgappen,
        maplim,minins,maxins,orient,rgid,rglb,rgpl,rgpu,rgsm,rgcn,rgds,rgdt,
        rgpi,rgpg,rgpm,dup,dupscore,dupregex,duppix,maxint,minreads,mismatches,
        windows,known,model,lod,entropy,maxcon,maxmoves,greedy,maxreads,
        maxinmem,covariates,knownsites,ics,maxcyc,mcs,bqsrpen,ddq,idq,calconf,
        emitconf,gtmode):
    project = os.path.dirname(os.path.abspath(__file__))
    ref_path = '/data/db/Homo_sapiens/UCSC/{0}/Sequence/Bowtie2Index/genome.fa'.format(reference)
    dbsnp = '/data/db/Homo_sapiens/UCSC/{0}/Annotation/Variation/dbsnp.vcf'.format(reference)
    config = configparser.ConfigParser()
    config['General'] = {'fwd':','.join(fwd), 'rev':','.join(rev), 'outdir': outdir,
                        'threads':threads, 'java':java, 'mem':mem, 
                        'reference': ref_path, 'dbsnp': dbsnp, 'pipeline':'dnaseq_cohort'}
    fastqc_path = '/projects/home/sravishankar9/tools/FastQC/fastqc'
    config['FastQC'] = {'fastqc': fastqc_path, 'kmer':kmer}
    trimmomatic = '{0}/Trimmomatic-0.33/trimmomatic-0.33.jar'.format(project)
    config['Trimmomatic'] = {'trimmomatic': trimmomatic, 'illuminaclip' : adap, 
                            'window': window, 'leading': lead, 'trailing': trail, 
                            'crop': crop, 'headcrop' : headcrop, 'minlen' : minlen}
    bowtiepath = '{0}/bowtie2-2.2.5/bowtie2'.format(project)
    bowtieindex = '/data/db/Homo_sapiens/UCSC/{0}/Sequence/Bowtie2Index/genome'.format(reference)
    config['Bowtie'] = {'bowtie': bowtiepath, 'index': bowtieindex, 'mode': '--'+mode, 
                        'mismatch': mismatch, 'seedlen':seedlen, 'interval': interval,
                        'seedex': seedex, 'reseed': reseed, 'ambfunc': ambfunc, 'dpad':dpad,
                        'gbar': gbar, 'matbonus': matbonus, 'penrange': penrange, 
                        'ambpen': ambpen, 'gappen': gappen, 'refgappen': refgappen, 
                        'maplim': maplim, 'minins': minins, 'maxins': maxins, 'orient': '--'+orient}
    picard = '{0}/picard-tools-1.139/picard.jar'.format(project)
    config['Picard'] = {'picard': picard, 'rgid':rgid, 'rglb': rglb, 'rgpl': rgpl, 
                        'rgpu' : rgpu, 'rgsm': rgsm, 'rgcn': rgcn, 'rgds': rgds, 'rgdt':rgdt,
                        'rgpi' : rgpi, 'rgpg' : rgpg, 'rgpm': rgpm, 'dup': dup, 
                        'dupscore': dupscore, 'dupregex': dupregex, 'duppix': duppix}
    gatk = '{0}/GenomeAnalysisTK/GenomeAnalysisTK.jar'.format(project)
    config['GATK'] = {'gatk': gatk, 'maxins': maxint, 'minreads': minreads,
                      'mismatch': mismatches, 'window': windows, 'known' : known,
                      'model': model, 'lod': lod, 'entropy': entropy, 'maxcon': maxcon,
                      'maxmoves': maxmoves, 'greedy': greedy, 'maxreads': maxreads, 
                      'maxinmem': maxinmem, 'covariates': covariates, 'knownsites': knownsites,
                      'ics': ics, 'maxcyc': maxcyc, 'mcs': mcs, 'bqsrpen': bqsrpen, 
                      'ddq': ddq, 'idq': idq, 'calconf': calconf, 'emitconf': emitconf,
                      'gtmode': gtmode}
    with open('config.cfg','w') as configfile:
        config.write(configfile)
    return(os.path.abspath(project+'/config.cfg'))

if __name__ == '__main__' :
    nester = argparse.ArgumentParser()
    nester.add_argument('-i', '--indir', type=str, help='Input directory path')
    nester.add_argument('-c', '--config', type=str, default=os.path.dirname(os.path.abspath(__file__))+'.config.cfg', help='Config file path')
    nester.add_argument('--outdir', type=str, default=None, help='Output directory path')
    nester.add_argument('--thread', type=str, default='6', help='Number of threads')
    nester.add_argument('--java', type=str, default='/usr/bin/java', help='Java path')
    nester.add_argument('--mem', type=str, default='4g', help='java memory limit')
    nester.add_argument('--reference', type=str, default='hg19', help='Genome build to use for analysis')
    nester.add_argument('--fastqckmer', type=str, default='7', help='FastQC kmer threshold')
    nester.add_argument('--window', type=str, default='4:15', help='Trimmomatic sliding window and quality threshold')
    nester.add_argument('--cropl', type=str, default='0', help='Crop 5\'')
    nester.add_argument('--cropt', type=str, default='0', help='Crop 3\'')
    nester.add_argument('--headcrop', type=str, default='0', help='Clip n bases from the 5\' end')
    nester.add_argument('--minlen', type=str, default='36', help='Minimum read length to retain read')
    nester.add_argument('--adapters', type=str, default='adapters.txt', help='FastQC adapter file')
    nester.add_argument('--mode', type=str, default='end-to-end', help='Bowtie2 alignment mode')
    nester.add_argument('--mismatch', type=str, default='0', help='Mismatches allowed')
    nester.add_argument('--seedlen', type=str, default='20', help='Seed length')
    nester.add_argument('--interval', type=str, default='S,1,0.50', help='Seed extension function')
    nester.add_argument('--seedx', type=str, default='20', help='Seed extension')
    nester.add_argument('--reseed', type=str, default='3', help='Reseed value')
    nester.add_argument('--ambfunc', type=str, default='L,0,0.15', help='ambfunc')
    nester.add_argument('--dpad', type=str, default='15', help='dpad')
    nester.add_argument('--gbar', type=str, default='4', help='gbar')
    nester.add_argument('--matbonus', type=str, default='2', help='matbonus')
    nester.add_argument('--penrange', type=str, default='6,2', help='penrange')
    nester.add_argument('--ambpen', type=str, default='1', help='ambpen')
    nester.add_argument('--gappen', type=str, default='5,3', help='gappen')
    nester.add_argument('--refgappen', type=str, default='5,3', help='refgappen')
    nester.add_argument('--maplim', type=str, default='10', help='maplim')
    nester.add_argument('--minins', type=str, default='0', help='minins')
    nester.add_argument('--maxins', type=str, default='500', help='maxins')
    nester.add_argument('--orient', type=str, default='fr', help='orientation')
    nester.add_argument('--rgid', type=str, default='null', help='rgid')
    nester.add_argument('--rglb', type=str, default='dnaseq_cohort', help='rglb')
    nester.add_argument('--rgpl', type=str, default='Illumina', help='rgpl')
    nester.add_argument('--rgpu', type=str, default='HiSeq2500', help='rgpu')
    nester.add_argument('--rgsm', type=str, default='null', help='rgsm')
    nester.add_argument('--rgcn', type=str, default='VannbergLab', help='rgcn')
    nester.add_argument('--rgds', type=str, default='null', help='rgds')
    nester.add_argument('--rgdt', type=str, default='null', help='rgdt')
    nester.add_argument('--rgpi', type=str, default='null', help='rgpi')
    nester.add_argument('--rgpg', type=str, default='null', help='rgpg')
    nester.add_argument('--rgpm', type=str, default='null', help='rgpm')
    nester.add_argument('--dup', type=str, default='false', help='dup')
    nester.add_argument('--dupscore', type=str, default='SUM_OF_BASE_QUALITIES', help='dupscore')
    nester.add_argument('--dupregex', type=str, default='[a-zA-Z0-9]+:[0-9]:([0-9]+):(0-9]+).*', help='dupregex')
    nester.add_argument('--duppix', type=str, default='100', help='duppix')
    nester.add_argument('--maxint', type=str, default='500', help='maxint')
    nester.add_argument('--minreads', type=str, default='4', help='minreads')
    nester.add_argument('--mismatches', type=str, default='0.0', help='mismatches')
    nester.add_argument('--windows', type=str, default='10', help='windows')
    nester.add_argument('--known', type=str, default='null', help='known indels')
    nester.add_argument('--model', type=str, default='USE_READS', help='model')
    nester.add_argument('--lod', type=str, default='5.0', help='lod')
    nester.add_argument('--entropy', type=str, default='0.15', help='entropy')
    nester.add_argument('--maxcon', type=str, default='30', help='maxcon')
    nester.add_argument('--maxmoves', type=str, default='200', help='maxmoves')
    nester.add_argument('--greedy', type=str, default='120', help='greedy')
    nester.add_argument('--maxreads', type=str, default='20000', help='maxreads')
    nester.add_argument('--maxinmem', type=str, default='150000', help='maxinmem')
    nester.add_argument('--covariantes', type=str, nargs='+', 
                        default='ReadGroupCovariate,QualityScoreCovariate,CycleCovariate,ContextCovariate', help='Covariates')
    nester.add_argument('--knownSites', type=str, nargs='+', default='None', help='Known sites')
    nester.add_argument('--ics', type=str, default='3', help='ics')
    nester.add_argument('--maxcyc', type=str, default='500', help='maxcyc')
    nester.add_argument('--mcs', type=str, default='2', help='mcs')
    nester.add_argument('--bqsrgpen', type=str, default='40.0', help='bqsrgen')
    nester.add_argument('--ddq', type=str, default='45', help='ddq')
    nester.add_argument('--idq', type=str, default='45', help='idq')
    nester.add_argument('--calconf', type=str, default='30.0', help='call confidence')
    nester.add_argument('--emitconf', type=str, default='30.0', help='emit confidence')
    nester.add_argument('--gtmode', type=str, default='DISCOVERY', help='genotype')
    
    args = nester.parse_args()
    if args.outdir == None:
        args.outdir = args.indir + 'Outputs'
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    read = re.compile('_[Rr]*1')
    basename = re.compile('_[Rr]*[12]')
    fastq = re.compile('fastq|fq')
    input_files = glob.glob(os.path.abspath(args.indir+'*.fastq'))
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
            fwd.append(input_dict[lines][0])
            rev.append(input_dict[lines][1])
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
            args.minreads, args.mismatches, args.windows, args.known, 
            args.model, args.lod, args.entropy, args.maxcon, args.maxmoves,
            args.greedy, args.maxreads,args.maxinmem, args.covariantes,
            args.knownSites, args.ics, args.maxcyc, args.mcs, args.bqsrgpen,args.ddq,
            args.idq, args.calconf,args.emitconf,args.gtmode)
    main()
