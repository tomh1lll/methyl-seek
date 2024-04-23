#####################################################################################################
# Bisulphite sequencing (methylseq) analysis workflow
#
# This Snakefile is designed to perform QC of raw bisulphite sequencing data,
# align reads to human (hg38) bisulfite genome using either Bismark (bowtie2) or
# bwa_meth, summarizes alignment stats, extracts methylated loci from the genome.
#
# Created: December 21, 2021
# Contact details: Tom Hill (tom.hill@nih.gov), Neelam Redekar (neelam.redekar@nih.gov),
# Skyler Kuhn (skyler.kuhn@nih.gov), Asya Khleborodova (asya.khleborodova@nih.gov)
#
# Last Modified: June 2, 2022
#
#####################################################################################################

from os.path import join
from snakemake.io import expand, glob_wildcards
from snakemake.utils import R
from os import listdir
import pandas as pd

##
## Locations of working directories and reference genomes for analysis
##
sample_file= config["samples"]
#rawdata_dir= config["rawdata_dir"]
working_dir= config["result_dir"]
hg38_fa= config["hg38_fa"]
phage_fa= config["phage_fa"]
hg38_gtf= config["hg38_gtf"]
hg38_rRNA_intervals= config["hg38_rRNA_intervals"]
hg38_bed_ref= config["hg38_bed_ref"]
hg38_refFlat= config["hg38_refFlat"]
bisulphite_genome_path= config["bisulphite_ref"]
phage_genome_path= config["phage_ref"]
bisulphite_fa= config["bisulphite_fa"]
species= config["species"]

REF_ATLAS=config["REF_ATLAS"]
CpG_MAP_TABLE=config["CpG_MAP_TABLE"]
REF_MARKERS=config["REF_MARKERS"]
REF_IDS=config["REF_IDS"]
GOLD_MARKERS=config["GOLD_MARKERS"]
HG38TOHG19=config["HG38TOHG19"]

## The file requires these headings as they are used in multiple rules later on.

## Here we read in the samples file generated and begin processing the data, printing out the samples and group comparisons

df = pd.read_csv(sample_file, header=None)
df1=df[0].str.split("/",expand=True)

SAMPLES=list(set(df1.iloc[: , -1].tolist()))
rawdata_dir="/".join(df1.iloc[0 , :-1].tolist()) 

#SAMPLES=list(set(df[0].tolist()))
#GROUPS=list(set(df['comp'].tolist()))

print(SAMPLES)
print(len(SAMPLES))
#print(GROUPS)
#print(len(GROUPS))

CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chr23','chrX']

RN = ['R1', 'R2']

rule All:
    input:
      # Creating data links:
      expand(join(working_dir, "raw/{samples}.{rn}.fastq.gz"), samples=SAMPLES, rn=RN),

      # Checking data quality:
      expand(join(working_dir, "rawQC/{samples}.{rn}_fastqc.html"), samples=SAMPLES, rn=RN),

      # bisulphite genome preparation
      join(bisulphite_genome_path, species, "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"),
      join(bisulphite_genome_path, species, "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"),

      # bismark align to human reference genomes
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.cram"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.flagstat"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.flagstat"),samples=SAMPLES),

      # get alignment statistics
      expand(join(working_dir, "bismarkAlign/{samples}.RnaSeqMetrics.txt"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.flagstat.concord.txt"),samples=SAMPLES),
      expand(join(working_dir, "rseqc/{samples}.inner_distance_freq.txt"),samples=SAMPLES),
      expand(join(working_dir, "rseqc/{samples}.strand.info"),samples=SAMPLES),
      expand(join(working_dir, "rseqc/{samples}.Rdist.info"),samples=SAMPLES),
      expand(join(working_dir, "QualiMap/{samples}/qualimapReport.html"),samples=SAMPLES),
      expand(join(working_dir, "QualiMap/{samples}/genome_results.txt"),samples=SAMPLES),
      expand(join(working_dir, "preseq/{samples}.ccurve"),samples=SAMPLES),
      expand(join(working_dir, "bismarkAlign/{samples}.star.duplic"), samples=SAMPLES),

      # extract CpG profile with bismark
      expand(join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),samples=SAMPLES),
      expand(join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.bismark_bt2_pe.deduplicated.bedGraph.gz"),samples=SAMPLES),

      # generate multiqc output
      "multiqc_report.html",

      # Deconvolution output
      expand(join(working_dir, "CpG_CSV/{samples}.csv"),samples=SAMPLES),
      expand(join(working_dir, "deconvolution_CSV/{samples}.csv"),samples=SAMPLES),
      expand(join(working_dir, "deconvolution_CSV/{samples}_deconv.log"),samples=SAMPLES),
      expand(join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.cfDNAmeInput.bedGraph"),samples=SAMPLES),
      expand(join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.cfDNAmeDeconvolution.tsv"),samples=SAMPLES),
      join(working_dir, "deconvolution_CSV/total.csv"),
      join(working_dir, "deconvolution_CSV/total_deconv_output.csv"),
      join(working_dir, "deconvolution_CSV/total_deconv_plot.png"),
      join(working_dir,"UXM/UXM_deconv.250.csv"),
      expand(join(working_dir,"UXM/{samples}.pat.gz"),samples=SAMPLES),
      expand(join(working_dir,"UXM/{samples}_deconv.250.csv"),samples=SAMPLES),


## Copy raw data to working directory
rule raw_data_links:
    input:
      join(rawdata_dir, "{samples}.{rn}.fastq.gz")
    output:
      join(working_dir, "raw/{samples}.{rn}.fastq.gz")
    params:
      rname="raw_data_links",
      dir=directory(join(working_dir, "raw")),
    shell:
      """
      mkdir -p {params.dir}
      ln -s {input} {output}
      """

## Run fastqc on raw data to visually assess quality
rule raw_fastqc:
    input:
      join(working_dir, "raw/{samples}.{rn}.fastq.gz")
    output:
      join(working_dir, "rawQC/{samples}.{rn}_fastqc.html")
    params:
      rname="raw_fastqc",
      dir=directory(join(working_dir, "rawQC")),
      batch='--cpus-per-task=2 --mem=8g --time=8:00:00',
    threads:
      2
    shell:
      """
      module load fastqc/0.11.9
      mkdir -p {params.dir}
      fastqc -o {params.dir} -f fastq --threads {threads} --extract {input}
      """


## trimming/filtering with trimGalore
rule trimGalore:
    input:
      F1=join(working_dir, "raw/{samples}.R1.fastq.gz"),
      F2=join(working_dir, "raw/{samples}.R2.fastq.gz"),
    output:
      temp(join(working_dir, "trimGalore/{samples}_val_1.fq.gz")),
      temp(join(working_dir, "trimGalore/{samples}_val_2.fq.gz")),
    params:
      rname="trimGalore",
      dir=directory(join(working_dir, "trimGalore")),
      tag='{samples}',
      fastqcdir=directory(join(working_dir, "postTrimQC")),
      command="--fastqc --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 --length 50 --gzip",
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=25g --time=10:00:00',
    threads:
      16
    shell:
      """
      module load trimgalore/0.6.7
      module load fastqc/0.11.9
      mkdir -p {params.dir}
      trim_galore --paired --cores {threads} {params.command} --basename {params.tag} --output_dir {params.dir} --fastqc_args "--outdir {params.fastqcdir}"  {input.F1} {input.F2}
      """

################## New edition - ended


## prepare bi-sulphite genome
rule prep_bisulphite_genome:
    input:
      bisulphite_fa
    output:
      join(bisulphite_genome_path, species, "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"),
      join(bisulphite_genome_path, species, "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"),
    params:
      rname="prep_bisulphite_genome",
      dir=directory(join(bisulphite_genome_path, species)),
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=20g --time=2:00:00',
    threads:
      16
    shell:
      """
      module load bismark/0.23.0
      mkdir -p {params.dir}
      cd {params.dir}
      #cp reference_genome {params.dir}/genome.fa
      bismark_genome_preparation --verbose --parallel {threads} {params.dir} #--single_fasta
      """

rule bismark_align:
    input:
      F1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
      F2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
    output:
      B1=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam")),
      B2=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.flagstat"),
      FQ1=temp(join(working_dir, "bismarkAlign/{samples}_val_1.fq.gz_unmapped_reads_1.fq.gz")),
      FQ2=temp(join(working_dir, "bismarkAlign/{samples}_val_2.fq.gz_unmapped_reads_2.fq.gz")),
      FQ3=temp(join(working_dir, "bismarkAlign/{samples}_val_1.fq.gz_ambiguous_reads_1.fq.gz")),
      FQ4=temp(join(working_dir, "bismarkAlign/{samples}_val_2.fq.gz_ambiguous_reads_2.fq.gz")),
    params:
      rname="bismark_align",
      dir=directory(join(working_dir, "bismarkAlign")),
      genome_dir=directory(join(bisulphite_genome_path, species)),
      command="--bowtie2 -N 1 --bam -L 22 --X 1000 --un --ambiguous -p 4 --score_min L,-0.6,-0.6",
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=100g --time=10:00:00',
      outbam=join(working_dir, "bismarkAlign/{samples}_val_1_bismark_bt2_pe.bam"),
      R1=join(working_dir,"bismarkAlign/{samples}_val_1_bismark_bt2_PE_report.txt"),
      R2=join(working_dir,"bismarkAlign/{samples}.bismark_bt2_PE_report.txt"),
    threads:
      16
    shell:
      """
      module load bismark/0.23.0 samtools
      mkdir -p {params.dir}
      bismark --multicore {threads} --temp_dir /lscratch/$SLURM_JOBID/ {params.command} --output_dir {params.dir} --genome {params.genome_dir} -1 {input.F1} -2 {input.F2}
      mv {params.outbam} {output.B1}
      mv {params.R1} {params.R2}
      samtools flagstat -@ {threads} {output.B1} > {output.B2}
      """

rule bismark_dedup:
    input:
      F1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),
    output:
      T1=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.deduplicated.bam")),
      B1=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam")),
      B2=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.flagstat"),
      C1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.cram"),
    params:
      rname="bismark_dedup",
      dir=directory(join(working_dir, "bismarkAlign")),
      prefix=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe"),
      genome=bisulphite_fa,
    threads:
      16
    shell:
      """
      module load bismark/0.23.0
      module load samtools/1.15
      cd {params.dir}
      samtools view -hb {input.F1} | samtools sort -n -@ {threads} -O BAM -o {output.T1}
      deduplicate_bismark --paired --bam --outfile {params.prefix} {output.T1}
      samtools view -T {params.genome} -@ 16 -h -C {output.B1} > {output.C1}
      samtools flagstat -@ {threads} {output.B1} > {output.B2}
      """

rule bismark_extract:
  input:
    bam=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    cov=join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),
    bed=join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.bismark_bt2_pe.deduplicated.bedGraph.gz"),
    CHGOB=temp(join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/CHG_OB_{samples}.bismark_bt2_pe.deduplicated.txt.gz")),
    CHGOT=temp(join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/CHG_OT_{samples}.bismark_bt2_pe.deduplicated.txt.gz")),
    CHHOB=temp(join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/CHH_OB_{samples}.bismark_bt2_pe.deduplicated.txt.gz")),
    CHHOT=temp(join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/CHH_OT_{samples}.bismark_bt2_pe.deduplicated.txt.gz")),
    CPGOB=temp(join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/CpG_OB_{samples}.bismark_bt2_pe.deduplicated.txt.gz")),
    CPGOT=temp(join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/CpG_OT_{samples}.bismark_bt2_pe.deduplicated.txt.gz")),
  params:
    rname='pl:bismark_extract',
    bismark_index=join(bisulphite_genome_path,species),
    outdir=join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated"),
    genome=bisulphite_fa,
  shell:
    """
    mkdir -p {params.outdir}
    module load bismark samtools bowtie
    bismark_methylation_extractor --paired-end --no_overlap --multicore 8 --gzip --report --bedGraph --counts --buffer_size 100G --no_header --cytosine_report --output {params.outdir} --scaffolds --genome_folder {params.bismark_index} {input.bam}
    """

rule prep_bisulphite_phage_genome:
    input:
      phage_fa
    output:
      join(phage_genome_path, "Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"),
      join(phage_genome_path, "Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"),
    params:
      rname="prep_phage_genome",
      dir=directory(phage_genome_path),
      batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=20g --time=2:00:00',
    threads:
      16
    shell:
      """
      module load bismark/0.23.0
      mkdir -p {params.dir}
      cd {params.dir}
      cp reference_genome {params.dir}/genome.fa
      bismark_genome_preparation --verbose --parallel {threads} {params.dir} #--single_fasta
      """

rule bismark_phage:
  input:
    F1=join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),
    F2=join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),
  output:
    join(working_dir, "bismark_phage/{samples}_val_1_bismark_bt2_pe.bam"),
  params:
    rname="bismark_phage",
    dir=directory(join(working_dir, "bismark_phage")),
    genome_dir=directory(phage_genome_path),
    command="--bowtie2 -N 1 --bam -L 22 --X 1000 --un --ambiguous -p 2 --score_min L,-0.6,-0.6",
    batch='--cpus-per-task=16 --partition=norm --gres=lscratch:100 --mem=100g --time=10:00:00',
  threads:
    16
  shell:
    """
    module load bismark/0.23.0
    mkdir -p {params.dir}
    bismark --multicore {threads} --temp_dir /lscratch/$SLURM_JOBID/ {params.command} --output_dir {params.dir} --genome {params.genome_dir} -1 {input.F1} -2 {input.F2}
    """

################## New edition - start
rule picard:
  input:
    file1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    bam=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bam")),
    bai=temp(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bai")),
    metrics=join(working_dir, "bismarkAlign/{samples}.star.duplic")
  params:
    rname='pl:picard',
    sampleName="{samples}",
  threads: 6
  shell:
    """
    module load samtools/1.15
    module load picard/2.26.9

    java -Xmx110g -XX:ParallelGCThreads=5 -jar ${{PICARDJARPATH}}/picard.jar AddOrReplaceReadGroups \
    I={input.file1} O=/lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.bam \
    TMP_DIR=/lscratch/$SLURM_JOBID RGID=id VALIDATION_STRINGENCY=SILENT RGLB=library RGPL=illumina RGPU=machine RGSM=sample;

    java -Xmx110g -XX:ParallelGCThreads=5 -jar ${{PICARDJARPATH}}/picard.jar MarkDuplicates \
    I=/lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.bam \
    O=/lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.dmark.bam \
    TMP_DIR=/lscratch/$SLURM_JOBID CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE={output.metrics};

    mv /lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.dmark.bam {output.bam};
    mv /lscratch/$SLURM_JOBID/{params.sampleName}.bismark_bt2_pe.dedup_rg_added.dmark.bai {output.bai};
    sed -i 's/MarkDuplicates/picard.sam.MarkDuplicates/g' {output.metrics};

    """

rule preseq:
  input:
    bam=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bam"),
  output:
    ccurve = join(working_dir, "preseq/{samples}.ccurve"),
  params:
    rname = "pl:preseq",
    dir = directory(join(working_dir, "preseq")),
  shell:
    """
    module load preseq/3.1.2
    mkdir -p {params.dir}
    preseq c_curve -B -o {output.ccurve} {input.bam}
    """

rule qualibam:
  input:
    bamfile=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.dedup_rg_added.dmark.bam"),
  output:
    report=join(working_dir,"QualiMap/{samples}/qualimapReport.html"),
    results=join(working_dir,"QualiMap/{samples}/genome_results.txt"),
  params:
    rname='pl:qualibam',
    outdir=join(working_dir,"QualiMap/{samples}"),
    gtfFile=hg38_gtf,
  threads: 16
  shell:
    """
    module load qualimap/2.2.1
    mkdir -p {params.outdir}
    unset DISPLAY;
    qualimap bamqc -bam {input.bamfile} --feature-file {params.gtfFile} -outdir {params.outdir} -nt {threads} --java-mem-size=120G
    """

rule rseqc:
  input:
    file1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    out1=join(working_dir,"rseqc/{samples}.strand.info"),
    out2=join(working_dir,"rseqc/{samples}.Rdist.info")
  params:
    bedref=hg38_bed_ref,
    dir=join(working_dir,"rseqc"),
    rname="pl:rseqc",
  shell:
    """
    module load rseqc/4.0.0
    mkdir -p {params.dir}
    infer_experiment.py -r {params.bedref} -i {input.file1} -s 1000000 > {output.out1}
    read_distribution.py -i {input.file1} -r {params.bedref} > {output.out2}
    """

rule inner_distance:
  input:
    bam=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    innerdists=join(working_dir,"rseqc/{samples}.inner_distance_freq.txt"),
  params:
    prefix=join(working_dir,"rseqc/{samples}"),
    dir=join(working_dir,"rseqc"),
    genemodel=hg38_bed_ref,
    rname="pl:inner_distance",
  shell:
    """
    module load rseqc/4.0.0
    mkdir -p {params.dir}
    inner_distance.py -i {input.bam} -r {params.genemodel} -k 10000000 -o {params.prefix}
    """

rule stats:
  input:
    file1=join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    outstar1=join(working_dir,"bismarkAlign/{samples}.RnaSeqMetrics.txt"),
    outstar2=join(working_dir,"bismarkAlign/{samples}.flagstat.concord.txt"),
  params:
    rname='pl:stats',
    refflat=hg38_refFlat,
    rrnalist=hg38_rRNA_intervals,
    picardstrand="SECOND_READ_TRANSCRIPTION_STRAND",
    statscript=join("workflow", "scripts", "bam_count_concord_stats.py"),
  shell:
    """
    module load python/3.8 samtools/1.15 picard/2.26.9
    java -Xmx110g -jar ${{PICARDJARPATH}}/picard.jar CollectRnaSeqMetrics REF_FLAT={params.refflat} I={input.file1} O={output.outstar1} RIBOSOMAL_INTERVALS={params.rrnalist} STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND TMP_DIR=/lscratch/$SLURM_JOBID  VALIDATION_STRINGENCY=SILENT;
    sed -i 's/CollectRnaSeqMetrics/picard.analysis.CollectRnaSeqMetrics/g' {output.outstar1}
    samtools flagstat {input.file1} > {output.outstar2};
    ## python3 {params.statscript} {input.file1} >> {output.outstar2} ## does require "Rstat" not available in biowulf
    """

################## New edition - ended

rule multiqc:
  input:
    expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.bam"),samples=SAMPLES),
    expand(join(working_dir, "bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),samples=SAMPLES),
    expand(join(working_dir, "trimGalore/{samples}_val_1.fq.gz"),samples=SAMPLES),
    expand(join(working_dir, "trimGalore/{samples}_val_2.fq.gz"),samples=SAMPLES),
  output:
    "multiqc_report.html",
  params:
    rname='pl:multiqc',
    dir=working_dir,
    bis_dir=directory(join(working_dir,"bismarkAlign")),
    script_dir=join(working_dir,"scripts"),
  shell:
    """
    module load multiqc/1.9 bismark
    cd {params.bis_dir}
    bismark2report
    bismark2summary
    cd {params.dir}
    multiqc --ignore '*/.singularity/*' -f --interactive .
    """

############### Deconvolution rules begin here
rule get_CpG:
  input:
    join(working_dir, "CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.bismark_bt2_pe.deduplicated.CpG_report.txt.gz"),
  output:
    join(working_dir, "CpG_CSV/{samples}.csv"),
  params:
    rname="get_CpG",
    cutoff=5,
    script_dir=join(working_dir,"scripts"),
    dir1=join(working_dir,"CpG_CSV"),
    dir2=join(working_dir,"deconvolution_CSV"),
  shell:
    """
    mkdir -p {params.dir1}
    mkdir -p {params.dir2}
    module load R
    Rscript {params.script_dir}/get_methy.R {input} {wildcards.samples} {params.cutoff} {output}
    """

rule get_overlap_meth:
  input:
    join(working_dir, "CpG_CSV/{samples}.csv"),
  output:
    join(working_dir, "deconvolution_CSV/{samples}.csv"),
  params:
    rname="get_overlap_meth",
    map_table=CpG_MAP_TABLE,
  run:
    df_ref=pd.read_csv(params.map_table,sep='\t',header=None)
    df_ref.columns=['chromosome','start','end','cgid']
    df_ref=df_ref.loc[(df_ref['chromosome'].isin(CHRS)),]
    dfm=pd.read_csv(input[0])
    dfm=pd.merge(df_ref,dfm,on=['chromosome','start','end'],how='inner')
    dfm=dfm.drop(labels=['chromosome','start','end'],axis=1)
    dfm=dfm.set_index('cgid')
    dfm.to_csv(output[0])

rule run_deconv:
  input:
    join(working_dir, "deconvolution_CSV/{samples}.csv"),
  output:
    join(working_dir, "deconvolution_CSV/{samples}_deconv.log"),
  params:
    script_dir=join(working_dir,"scripts"),
    dir=join(working_dir,"deconvolution_CSV"),
    rname="run_deconv",
    ref=REF_ATLAS,
  shell:
    """
    module load python
    cd {params.dir}
    python {params.script_dir}/deconvolve.py --atlas_path {params.ref} --plot --residuals {input}  > {output}  2>&1
    """

rule merge_tables:
  input:
    expand(join(working_dir, "deconvolution_CSV/{samples}.csv"),samples=SAMPLES),
  output:
    join(working_dir, "deconvolution_CSV/total.csv"),
  params:
    rname="merge_tables",
  run:
    dfm=pd.read_csv(input[0])
    for f in input[1:]:
        df=pd.read_csv(f)
        dfm=pd.merge(dfm,df,on='cgid',how='outer')
    dfm.to_csv(output[0],index=False)

rule run_deconv_merged:
  input:
    join(working_dir, "deconvolution_CSV/total.csv"),
  output:
    join(working_dir, "deconvolution_CSV/total_deconv_output.csv"),
    join(working_dir, "deconvolution_CSV/total_deconv_plot.png"),
  params:
    ref=REF_ATLAS,
    dir=join(working_dir,"deconvolution_CSV"),
    script_dir=join(working_dir,"scripts"),
    rname="run_deconv_merged",
  shell:
    """
    module load python
    cd {params.dir}
    python {params.script_dir}/deconvolve.py --atlas_path {params.ref} --plot --residuals {input}
    """



###############CFDNAME RULES######################

rule sorting_CpG_bedgraph:
  input:
    bed=join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.bismark_bt2_pe.deduplicated.bedGraph.gz"),
  output:
    bed=temp(join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.mapped_autosomal_CpG.bedGraph.tmp")),
  params:
    rname="pl:sorting_CpG_bedgraph",
  shell:
    """
    module load bedtools
    zcat {input.bed} | tail -n +2 | bedtools sort -i - > {output.bed}
    """

rule liftover_bedgraph:
  input:
    bed=join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.mapped_autosomal_CpG.bedGraph.tmp"),
  output:
    graph=temp(join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.liftover.bedGraph.tmp")),
  params:
    rname="pl:liftover_bedgraph",
    lift_file=HG38TOHG19,
  shell:
    """
    module load bedtools crossmap/0.6.3
    crossmap bed {params.lift_file} {input.bed} {output.graph}
    """

rule extract_signature_beds:
  input:
    graph=join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.liftover.bedGraph.tmp"),
  output:
    sort=temp(join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.sorted.bedGraph")),
  params:
    rname="pl:extract_signature_beds",
    markers=GOLD_MARKERS,
  shell:
    """
    module load bedtools
    bedtools sort -i {input.graph} | bedtools intersect -wo -a {params.markers} -b stdin -sorted | awk '$6-$5==1 {{print $0}}' | awk 'NF{{NF-=1}};1' > {output.sort}
    """

rule aggregate_over_regions:
  input:
    sort=join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.sorted.bedGraph"),
  output:
    tsv=join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.cfDNAmeInput.bedGraph"),
  params:
    rname="pl:aggregate_over_regions",
    script_dir=join(working_dir,"scripts"),
  shell:
    """
    module load R
    Rscript {params.script_dir}/aggregate_over_regions.R {input.sort} {output.tsv}
    """

rule cfDNAme:
  input:
    bed=join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.cfDNAmeInput.bedGraph"),
  output:
    tsv=join(working_dir,"CpG/{samples}.bismark_bt2_pe.deduplicated/{samples}.cfDNAmeDeconvolution.tsv"),
  params:
    rname="pl:cfDNAme",
    script_dir=join(working_dir,"scripts"),
    sampleName="{samples}",
    reference_markers=REF_MARKERS,
    reference_IDs=REF_IDS,
  shell:
    """
    module load R
    Rscript ${params.script_dir}/tissues_of_origin_v2.R \
    {input.bed}  \
    {params.reference_markers} \
    {output.tsv} \
    {params.reference_IDs} \
    {params.sampleName} \
    TRUE \
    FALSE \
    TRUE \
    colon_5/hsc_5/hsc_2/spleen_1/spleen_2/spleen_3/spleen_4/
    """

rule bamsort:
  input:
    bam=join(working_dir,"bismarkAlign/{samples}.bismark_bt2_pe.deduplicated.bam"),
  output:
    bam=temp(join(working_dir,"bismarkAlign/{samples}.bam")),
    bai=temp(join(working_dir,"bismarkAlign/{samples}.bam.bai")),
  params:
    rname="pl:bamsort",
    reference="/data/NHLBI_IDSS/references/Bismark_Genomes/hg38/genome.fa",
  shell:
    """
    module load samtools 
    samtools sort -@ 12 --reference /data/NHLBI_IDSS/references/Bismark_Genomes/hg38/genome.fa {input.bam} > {output.bam}
    samtools index -@ 12 {output.bam}
    """

rule wgbstools:
  input:
    bam=join(working_dir,"bismarkAlign/{samples}.bam"),
    bai=temp(join(working_dir,"bismarkAlign/{samples}.bam.bai")),
  output:
    pat=join(working_dir,"UXM/{samples}.pat.gz"),
  params:
    rname="pl:wgbstools",
    ref=species,
    outdir=join(working_dir,"UXM"),
  shell:
    """
    export PATH=${PATH}:/data/NHLBI_IDSS/references/UXM
    module load samtools bedtools bamtools
    wgbstools bam2pat --genome {params.ref} --out_dir {params.outdir} -@ 12 {input.bam}
    """

rule UXM:
  input:
    pat=join(working_dir,"UXM/{samples}.pat.gz"),
  output:
    pat=join(working_dir,"UXM/{samples}_deconv.250.csv"),
  params:
    rname="pl:UXM",
    atlas="/data/NHLBI_IDSS/references/UXM/supplemental/Atlas.U250.l4.hg38.full.tsv",
  shell:
    """
    export PATH=${{PATH}}:/data/NHLBI_IDSS/references/UXM
    module load samtools bedtools bamtools
    uxm deconv {input.pat} -o {output.pat} --atlas {params.atlas} --ignore Colon-Fibro Dermal-Fibro Gallbladder Bone-Osteob
    """

rule UXM_all:
  input:
    pat=expand(join(working_dir,"UXM/{samples}.pat.gz"),samples=SAMPLES),
  output:
    pat=join(working_dir,"UXM/UXM_deconv.250.csv"),
  params:
    rname="pl:UXM_all",
    atlas="/data/NHLBI_IDSS/references/UXM/supplemental/Atlas.U250.l4.hg38.full.tsv",
  shell:
    """
    export PATH=${{PATH}}:/data/NHLBI_IDSS/references/UXM
    module load samtools bedtools bamtools
    uxm deconv {input.pat} -o {output.pat} --atlas {params.atlas} --ignore Colon-Fibro Dermal-Fibro Gallbladder Bone-Osteob
    """
