##NEST

NEST is a pipelining toolkit allowing users to run multiple gold standard pipelines on NGS Data. It designed to be a self sufficient framework, requiring minimal environment setup. Using python's multiprocessing capabiities, it can be used to analyze multiple samples and handle biological as well technical replicates. The parameters can either be set using a configuration script or by directly editing the config file. The framework allows you to customize your analysis and also restart analysis from any point in the pipeline.

###Release Notes:
* DNA/Exome Seq analysis pipeline (GATK gold standard).
* Annotation using snpEff.
* Quality trimming using TrimGalore.
* Adapters and contaminant can be supplied as external files to ensure accurate read trimming.

###Software Requirements:
* [Python 3.4](https://www.python.org/download/releases/3.4.1/)
 



