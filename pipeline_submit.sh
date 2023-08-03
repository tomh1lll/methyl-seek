#! /bin/bash
set -euo pipefail

###################
#
# Launching shell script for NHLBI processing of methylation data
#
###################
module load python/3.7
module load snakemake/5.13.0

cd $SLURM_SUBMIT_DIR

R=$3
INPUT_DIRECTORY=$2

find "$INPUT_DIRECTORY" -iname '*R?.fastq.gz' | sed 's/.R[1,2].fastq.gz//g' | sort | uniq > $R/samples.txt

echo $R

echo "generate config.yaml"

cat << EOF > $R/config.yaml
samples: "${R}/samples.txt"
result_dir: "${R}"
hg38_fa: "/data/NHLBI_IDSS/references/Bismark_Genomes/hg38/genome.fa"
hg38_gtf: "/fdb/igenomes/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
hg38_rRNA_intervals: "/data/NHLBI_IDSS/references/Bismark_Genomes/hg38/RNA-Seek_build/hg38.rRNA_interval_list"
hg38_refFlat: "/data/NHLBI_IDSS/references/Bismark_Genomes/hg38/RNA-Seek_build/refFlat.txt"
hg38_bed_ref: "/data/NHLBI_IDSS/references/Bismark_Genomes/hg38/RNA-Seek_build/genes.ref.bed"
mm10_fa: "/fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
rn6_fa: "/fdb/igenomes/Rattus_norvegicus/UCSC/rn6/Sequence/WholeGenomeFasta/genome.fa"
bisulphite_ref: "/data/NHLBI_IDSS/references/Bismark_Genomes"
bisulphite_fa: "/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
phage_fa: "/data/NHLBI_IDSS/references/Bismark_Genomes/phage/genome.fa"
phage_ref: "/data/NHLBI_IDSS/references/Bismark_Genomes/phage"
species: "hg38"
REF_ATLAS: "/data/NHLBI_IDSS/references/MethylCpG_deconvolution/resources/reference_atlas.csv"
CpG_MAP_TABLE: "/data/NHLBI_IDSS/references/MethylCpG_deconvolution/resources/reference_atlas.hg38.unique.bed"
REF_MARKERS: "/data/NHLBI_IDSS/projects/NHLBI-4/cfDNAme_data/Atlas.MethylMatrix.markers_22cov.tsv"
REF_IDS: "/data/NHLBI_IDSS/projects/NHLBI-4/cfDNAme_data/reference_methylomes.20230803.tsv"
GOLD_MARKERS: "/data/NHLBI_IDSS/projects/NHLBI-4/cfDNAme_data/golden_markers.22.bed"
HG38TOHG19: "/data/NHLBI_IDSS/projects/NHLBI-4/cfDNAme_data/hg38ToHg19.over.chain.gz"

EOF

mkdir -p $R/snakejobs
mkdir -p $R/reports

##
## Test commandline arguments
##

if [ $1 != "npr" ] && [ $1 != "process" ] ; then
    echo " "
    echo "Invalid commandline option: $1"
    echo "Valid commandline options include: npr or process"
    echo " "
    exit
fi

# and be sure this directory and all subdirectories are writable
#chmod -fR g+rwx .
#chmod +rwx .

SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`

echo $SCRIPT
echo $SCRIPTPATH

##
## Run snakemake
##
echo "Run snakemake"

CLUSTER_OPTS="sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name={params.rname} -e snakejobs/slurm-%j_{params.rname}.out -o snakejobs/slurm-%j_{params.rname}.out"

if [ $1 == "npr" ]
then
    snakemake -npr --snakefile $R/Snakefile -j 1 --configfile $R/config.yaml
fi

if [ $1 == "process" ]
then
    snakemake --latency-wait 120  -s $R/Snakefile -d $R --printshellcmds --configfile $R/config.yaml --cluster-config $R/cluster.json --keep-going --restart-times 1 --cluster "$CLUSTER_OPTS" -j 500 --rerun-incomplete --stats $R/reports/snakemake.stats | tee -a $R/reports/snakemake.log
fi

