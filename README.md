# cutrun
Ting Lab Nextflow pipeline for processing CUT&amp;RUN, ChIP-seq, and ChIP-like data

**WIP**: This is a current work in progress. The main pipeline should work, however spike-in implementation is currently a work in progress. The `spike-in` directory currently contains all the scripts used to create the T4phage bowtie2 indices, which are required for the spike-in steps of the pipeline. Please direct all questions towards bradlem4@ccf.org. 

## Software Requirements
The following packages are loaded using `module load` in the pipeline. If your platform does not support environment modules, or doesn't have modules for the following software then it is suggested that you install them via conda. Additionally, the version numbers have been included below for the version of the software that was tested, however there's no reason to believe a different version of the same software wouldn't work. 
```
samtools == 1.14
bowtie2 == 2.3.4.1
FastQC == 0.11.9
```
Additionally, `macs2`, `picard`, and `deeptools` are required. I have included lines in the respective processes to install these on the fly using `conda`, however if you would like to build your own environments you may simply replace the `conda` directive argument with the path to the conda environment on your file system ([as seen here](https://www.nextflow.io/docs/latest/conda.html))

You will also need to download a reference genome and indices (hg38.fa) as well as the corresponding `bowtie2` indices. In a later version of the pipeline, a step to download these may be included.
