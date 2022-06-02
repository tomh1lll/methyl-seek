
#### Setting up the working environment

Before running the pipeline, certain packages are required to be installed within a custom conda environment.

```
module load python
source /data/$USER/conda/etc/profile.d/conda.sh
conda create --name meth
conda activate meth
mamba install -yc bioconda bwameth methyldackel
conda deactivate meth
```

#### Setting up the working files

Alter the config.yaml so the rawdata_dir is the absolute path of the directory containing all your fastqs.
Alter the result_dir so it is the absolute path of the working directory containing your snakemake pipeline, where results will be stored.
Alter samples in config.yaml to be the absolute path of your samples.txt. Check this is correct. The samples file should have a list of all sample IDs, one per line, as follows:

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
