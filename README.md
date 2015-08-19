##Vannberg lab variant caller

This is python pipeline which integrates various open source tools to analyze whole exome and whole genome seqeuncing data. The pipeline takes in single end or paired end fastq files, performs quality checks, read alignment, indel realignment, variant calling and variant annotation to give a fully annotated vcf file. This current version also annotates the variants with information from COSMIC. 
The pipeline needs, other than the fastq files, the reference genome, the run information of the fastq files and filter parameters(defaults filter parameters are used if not specified).

###Software Requirements:
* [Python 2.7](https://www.python.org/download/releases/2.7/)
* [Pysam](https://github.com/pysam-developers/pysam)
* [Pyvcf](https://github.com/jamescasbon/PyVCF)
* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [GATK](https://www.broadinstitute.org/gatk/)
* [Samtools](http://samtools.sourceforge.net/)
* [Picard](http://sourceforge.net/projects/picard/)
* [Fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Sickle](https://github.com/najoshi/sickle)
* [Annovar](http://www.openbioinformatics.org/annovar/)
 



