

#### Setting up the working files

During the submission process, you provide the script with the directories where rawdata is stored & where the output will be stored. First make sure that these two directories exist & that the raw data directory contains all the fastqs of samples that would you to process. The fastqs need to be paired-end and named in the following format.

```
S1.R1.fastq.gz
S1.R2.fastq.gz
S2.R1.fastq.gz
S2.R2.fastq.gz
S3.R1.fastq.gz
S3.R2.fastq.gz
S4.R1.fastq.gz
S4.R2.fastq.gz
S5.R1.fastq.gz
S5.R2.fastq.gz
S6.R1.fastq.gz
S6.R2.fastq.gz
```

where S1-6 are arbitrarily chosen example sample names. The samples matching this name scheme in the rawdata directory will generate a sample file with the following names in the output directory:

```
S1
S2
S3
S4
S5
S6
```

The pipeline performs quality control steps and maps the data using Bismark, before extracting CpG profiles using Bismark. Following this the pipeline uses the previously generated CpG to deconvolute the data and identify which tissues samples belong to based on methylation profiles.

#### Dry run of the pipeline

To perform a dry run a step of the pipeline, giving 3 arguments, if it is a dry-run (npr) or full run (process), the directory containing all the samples you want to process, and the output working directory where results will be stored, and submit:

```
sh pipeline_submit.sh npr /data/NHLBIcore/rawdata/methyl-seek /data/NHLBIcore/projects/methyl-seek
```

#### Actual run of the pipeline

Once everything seems to work, to perform a full run of a step of the pipeline, submit:

```
sbatch --partition=norm --gres=lscratch:500 --time=10-00:00:00 --mail-type=BEGIN,END,FAIL pipeline_submit.sh process /data/NHLBIcore/rawdata/methyl-seek /data/NHLBIcore/projects/methyl-seek
```
