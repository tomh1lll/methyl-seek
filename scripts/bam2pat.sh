#!/bin/bash

export PATH=${PATH}:/data/NHLBI_IDSS/references/UXM

module load samtools bedtools bamtools

in=$1

out=$2

outdir=$3

mkdir -p ${outdir}

samtools sort -@ 12 --reference /data/NHLBI_IDSS/references/Bismark_Genomes/hg38/genome.fa ${in} > ${out}

samtools index ${out}

wgbstools bam2pat --genome hg38 --out_dir ${outdir} -@ 12 ${out}
