##NEST

NEST is a pipelining toolkit allowing users to run multiple gold standard pipelines on NGS Data. It designed to be a self sufficient framework, requiring minimal environment setup. Using python's multiprocessing capabiities, it can be used to analyze multiple samples and handle biological as well technical replicates. The parameters can either be set using a configuration script or by directly editing the config file. Allows for three analysis modes for both Exome and DNAseq:

![pipeline](Nest.png)




###Release Notes:

####Version0.9:
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
 



