# Docker image for MiXCR performance analysis with _in silico_ generated data

The Docker image contains pre-installed software and scripts to generate _in silico_ RNA-Seq data and run performance analysis.


# Running analysis

Docker image is already built and deposited in Docker Hub. To run the full analysis just install [Docker](https://www.docker.com) and execute the following line (no any special downloads required, Docker will automatically download the image):
 ```
 docker run -v WORKING_PATH:/work -i milaboratory/rna-seq-paper:v1.0 /bin/bash -c "/opt/scripts/run-comparison.sh && /opt/scripts/run-false-positives.sh" > log 2>errLog
 ```
where `WORKING_PATH` should be replaced with some particular path at your host computer. Be sure that you have at least ~80Gb of RAM and ~1TB of disk space.

# Description 

Scripts [`run-comparison.sh`](run-comparison.sh) and [`run-false-positives.sh`](run-false-positives.sh) will execute the following steps:
 - generate a set of TRB VDJ recombinations with [repseqio](https://github.com/repseqio/repseqio) tool
 - generate `fastq` files with _in silico_ RNA-seq data with a) a portion of TRB records b) without any TCR records
 - run [MiXCR](http://github.com/milaboratory/mixcr) analysis in RNA-Seq mode
 - run [STAR](https://github.com/alexdobin/STAR) aligner to produce BAMs and then run [TRUST](http://www.nature.com/ng/journal/v49/n4/full/ng.3820.html?WT.feed_name=subjects_systems-biology)
 - compare results obtained by MiXCR and TRUST with the original set of CDR3 clones and draw comparison plots
 - calculate rate of false-extensions produced by MiXCR `extendAlignments` action
 - calculate rate of false-overlaps produced by MiXCR `assemblePartial` action
 
All intermediate data can be found in the `WORKING_PATH`. 

The following resulting files will be produced:
 - `SupplementaryFigure_in_silico.pdf` &mdash;  comparison of MiXCR and TRUST results
 - `falseExtensionsResults.txt` &mdash;  MiXCR false-positive extension rate estimations
 - `falseOverlapsResults.txt` &mdash;  MiXCR false-positive partial assembly rate estimations
     
For the details on available options see [`run-comparison.sh`](run-comparison.sh) script.
