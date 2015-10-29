##NEST

NEST is a pipelining toolkit allowing users to run multiple gold standard pipelines on NGS Data. It designed to be a self sufficient framework, requiring minimal environment setup. Using python's multiprocessing capabiities, it can be used to analyze multiple samples and handle biological as well technical replicates. The parameters can either be set using a configuration script or by directly editing the config file. Allows for three analysis modes for both Exome and DNAseq:

![pipeline](Nest.png)


Running Nest:
```python
./nest.py -h
usage: nest.py [-h] [-i INDIR] [-c CONFIG] [--outdir OUTDIR] [--thread THREAD]
               [--java JAVA] [--mem MEM] [--reference REFERENCE]
               [--pipeline {exomeseq_cohort,dnaseq_cohort,dnaseq,exomeseq,exomeseq_rapid}]
               [--exome EXOME] [--fastqckmer FASTQCKMER] [--window WINDOW]
               [--cropl CROPL] [--cropt CROPT] [--headcrop HEADCROP]
               [--minlen MINLEN] [--adapters ADAPTERS] [--mode {local,global}]
               [--mismatch MISMATCH] [--seedlen SEEDLEN] [--interval INTERVAL]
               [--seedx SEEDX] [--reseed RESEED] [--ambfunc AMBFUNC]
               [--dpad DPAD] [--gbar GBAR] [--matbonus MATBONUS]
               [--penrange PENRANGE] [--ambpen AMBPEN] [--gappen GAPPEN]
               [--refgappen REFGAPPEN] [--maplim MAPLIM] [--minins MININS]
               [--maxins MAXINS] [--orient ORIENT] [--rgid RGID] [--rglb RGLB]
               [--rgpl RGPL] [--rgpu RGPU] [--rgsm RGSM] [--rgcn RGCN]
               [--rgds RGDS] [--rgdt RGDT] [--rgpi RGPI] [--rgpg RGPG]
               [--rgpm RGPM] [--dup DUP] [--dupscore DUPSCORE]
               [--dupregex DUPREGEX] [--duppix DUPPIX] [--maxint MAXINT]
               [--minreads MINREADS] [--mismatches MISMATCHES]
               [--windows WINDOWS] [--model MODEL] [--lod LOD]
               [--entropy ENTROPY] [--maxcon MAXCON] [--maxmoves MAXMOVES]
               [--greedy GREEDY] [--maxreads MAXREADS] [--maxinmem MAXINMEM]
               [--covariantes COVARIANTES [COVARIANTES ...]] [--ics ICS]
               [--maxcyc MAXCYC] [--mcs MCS] [--bqsrgpen BQSRGPEN] [--ddq DDQ]
               [--idq IDQ] [--calconf CALCONF] [--emitconf EMITCONF]
               [--gtmode GTMODE] [--mmq MMQ] [--erc ERC] [--ann ANN]
               [--applymode APPLYMODE] [--filterlevel FILTERLEVEL]
```
By default NEST looks for config file within the NEST directory.
Outputs folder stored in folder containing input files by default.

SeqMultiple DNASeq mode:
```python
./nest.py -i <input directory path> --pipeline dnaseq_multiple
```

SeqCohort DNASeq mode:
```python
./nest.py -i <input directory path> --pipeline dnaseq_cohort
```

SeqRapid DNASeq mode:
```python
./nest.py -i <input directory path> --pipeline dnaseq_rapid
```
SeqMultiple ExomeSeq mode:
```python
./nest.py -i <input directory path> --pipeline exomeseq_multiple
```

SeqCohort ExomeSeq mode:
```python
./nest.py -i <input directory path> --pipeline exomeseq_cohort
```

SeqRapid ExomeSeq mode:
```python
./nest.py -i <input directory path> --pipeline exomeseq_rapid
```


###Release Notes:

####Version0.9.1:
* SnpEFF replaced by Annovar.
* Mendelian case control analysis scripts added.
* Vcf plotting tools added.
* Illumina Rapid run feature added.
* VQSR added.
* Issue with memory and thread allocation per sample resolved.

####Version0.9.0:
* Introducing config files to track runs.
* Default config file generated each run, can be altered using command line options.
* User can supply a custom config file, to run analysis with desired settings.
* Single and multiple sample DNASeq/ ExomeSeq analysis.
* Dependencies packaged in the repository.
* Uses Illumina igenomes reference.
* Trimgalore replaced by Trimmomatic.
* UnifiedGenotyper replaced by HaplotypeCaller.

####Version 0.8.9:
* DNA/Exome Seq analysis pipeline (GATK gold standard).
* Annotation using snpEff.
* Quality trimming using TrimGalore.
* Adapters and contaminant can be supplied as external files to ensure accurate read trimming.

###Software Requirements:
* [Python 3.4](https://www.python.org/download/releases/3.4.1/)
* [Pandas](http://pandas.pydata.org/)
* [Seaborn](http://web.stanford.edu/~mwaskom/software/seaborn/index.html) 



