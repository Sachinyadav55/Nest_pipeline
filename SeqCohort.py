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

class SeqCohort:
    
    #Regex and dictionary initialization
    #For output directory sturcture creation

    def __init__(self):
        read = re.compile('_[Rr]*1')
        output_record = defaultdict()
        file_record = defaultdict()
        config = configparser.ConfigParser()
        config.read('config.cfg')
        #Get samples, sample names, reference, and outdir
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
        #Create sample specific output paths
        self.outpaths = [self.outdir + '/' + val for val in \
                        self.samples]
        output_record = dict()
        file_record = dict()
        for paths,fwd,rev,name in zip(self.outpaths,self.read1, self.read2,
            self.samples):
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
        self.fwd, self.rev = trimgalore.run_qc(self.fwd, self.rev, qcdir, 
            self.base)

        #Run alignment
        self.bamfile = bowtie.bowtie(self.fwd, self.rev, self.datadir, 
            self.base)
        
        #Run picard tools 
        self.bamfile = picard.addreadgroup(self.bamfile, self.base)
        self.bamfile, self.metrics = picard.markdup(self.bamfile)

        #Run GATK realignment and quality recalibration
        self.interval = gatk.target(self.bamfile)
        self.bamfile = gatk.realigner(self.bamfile, self.interval)
        self.table = gatk.baserecal(self.bamfile)
        self.bamfile = gatk.printreads(self.bamfile, self.table)
        self.gvcf = gatk.haplocaller(self.bamfile)
        return(self.gvcf)

    def joint_genotyping(self, gvcf_list):
        #Run joint genotyping
        self.gvcf_list = gvcf_list
        self.vcffile = gatk.joint_genotyper(gvcf_list, self.outdir)
        return(self.vcffile)

    def vqsr(self, vcffile):
        #Run VQSR and apply reclibration
        self.vcffile = vcffile
        self.srf, self.irf, self.stf, self.itf = gatk.vqsr(vcffile, self.outdir)
        self.vcffile = gatk.applyrecalibration(self.srf, self.irf, self.stf,
            self.itf, self.vcffile, self.outdir)
        return(self.vcffile)


