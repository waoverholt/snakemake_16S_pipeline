# 16S Amplicon Initial Processing Pipeline

This is a snakemake pipeline for initial handling of 16S amplicon datasets. It seeks to automate many of the steps that I describe [here](http://waoverholt.github.io/blog/).

It takes as raw input sequences in the format SAMPLE-NAME_R1_001.fastq.gz that our MiSeq produces by default. The final output is an unfiltered OTU table suitable for further work with R, MOTHUR, or QIIME. This repository will not likely be maintained and I will probably end up testing QIIME2 and might migrate to newer OTU picking methods. However, feel free to use it as a reference (crude) snakemake pipeline. Also, don't hesitate to contact me at waoverholt@gmail.com if you have specific questions about this pipeline.

# Dependencies
This pipeline was written with [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and requires python3. It is set up to be configurable with miniconda.
Follow instructions to install [python](https://www.python.org/downloads/) and [miniconda](https://conda.io/docs/user-guide/install/index.html).

Once you have miniconda installed you can use the environment.yaml file to install all dependencies in a virtual environment.
`conda env create -n snakemake_16S python=3.5 --file environment.yaml`
The environment name can be anything you'd like.

To set up your own environment, the pipeline requires the following packages and programs in your path:
 - [pear](https://sco.h-its.org/exelixis/web/software/pear/doc.html)
 - [vsearch](https://github.com/torognes/vsearch)
 - [mothur](https://www.mothur.org/)
 - [swarm](https://github.com/torognes/swarm)
 - [biom-format](http://biom-format.org/index.html)
 - [numpy](http://www.numpy.org/)
 - [scipy](https://www.scipy.org/)
 - [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
 - [pyyaml](https://github.com/yaml)

# Config File
The pipeline requires a config.yaml file to run. Please modify the existing config file for your datasets.

# Pipeline Summary
1. Merge paired reads with pear.
2. Quality control with vsearch.
3. Adapter triming with MOTHUR.
4. Dereplication with vsearch.
5. Denovo and reference chimera detection with vsearch.
6. OTU picking with swarm.

#Chimera reference database
The chimera reference database I use it too large to be hosted on github. I am currently using the 97 rep_set from [SILVA128 database](https://www.arb-silva.de/download/archive/qiime/).

The resulting OTU table is not abundance filtered. The tab delimited and the biom format tables produced are identical. Further analyzes can proceed using QIIME or R. To work with QIIME, you will need to deactivate the conda environment.
