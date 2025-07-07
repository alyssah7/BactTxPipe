#Snakefile:
import sys
import os

RIBOSOME_RNA_fasta = "VE303-04_rRNA_seqs_PROKKA.fasta"
ANNOT_REF_GENOME="/orcd/data/tami/003/projects/aro_vedanta/prokka/annotation/VE303-04_hybrid-v3/prokka_out.ffn"

PROJECT_SCRIPTS_DIRECTORY='./scripts'
sys.path.insert(0, PROJECT_SCRIPTS_DIRECTORY)
CURRENT_DIRECTORY = os.getcwd()
from gus_helper_functions import *
spls = 'samples.csv'

[PATH_ls, SAMPLE_ls, FILENAME_ls, REF_Genome_ls, GROUP_ls, OUTGROUP_ls] = read_samples_CSV(spls)
# Write sample_info.csv for each sample
split_samplesCSV(PATH_ls, SAMPLE_ls, FILENAME_ls, REF_Genome_ls, GROUP_ls, OUTGROUP_ls)
UNIQ_GROUP_ls = set(GROUP_ls)


rule all:
	input:
		# # Through lab mapping pipeline # #
		"1-bowtie2/alignment_stats.csv",
		"3-DESeq2/deseq_results.csv",
		

# link to sequence files
rule make_data_links:
  # NOTE: All raw data needs to be names fastq.gz. No fq! The links will be names fq though.
    input:
        sample_info_csv="data/{sampleID}/sample_info.csv",
    params:
        links_dir='links',
    output:
      # Recommend using symbolic links to your likely many different input files
        fq1="links/{sampleID}/R1.fq.gz",
        fq2="links/{sampleID}/R2.fq.gz",
    run:
        # get stuff out of mini csv file
        paths, sam, ref, fn = read_sample_info_CSV(input.sample_info_csv)
        print(paths)
        print(f"{sam},{ref},{fn}")
        # make links
        subprocess.run('mkdir -p links', shell=True)
        if len(paths)>1:
            cp_append_files(paths, sam, fn, params.links_dir)
        else:
            makelink(paths[0], sam, fn, params.links_dir)

### process reads ###
rule cutadapt:
  input:
      fq1 = "links/{sampleID}/R1.fq.gz",
      fq2 = "links/{sampleID}/R2.fq.gz",
  output:
      fq1o="tmp/{sampleID}_R1_trim.fq.gz",
      fq2o="tmp/{sampleID}_R2_trim.fq.gz",
  log:
      log="logs/cutadapt_{sampleID}.txt",
  conda:
      "cutadapt",
  shell:
      "cutadapt -a CTGTCTCTTAT --cores=4 -o {output.fq1o} {input.fq1} 1> {log};"
      "cutadapt -a CTGTCTCTTAT --cores=4 -o {output.fq2o} {input.fq2} 1>> {log};"

rule sickle2050:
  input:
      fq1o="tmp/{sampleID}_R1_trim.fq.gz",
      fq2o="tmp/{sampleID}_R2_trim.fq.gz",
  output:
      fq1o="tmp/{sampleID}_filt/filt1.fq.gz",
      fq2o="tmp/{sampleID}_filt/filt2.fq.gz",
      fqSo="tmp/{sampleID}_filt/filt_sgls.fq.gz",
  log:
      log="logs/sickle2050_{sampleID}.txt",
  conda:
      "envs/sickle-trim.yaml",
  shell:
      "sickle pe -g -f {input.fq1o} -r {input.fq2o} -t sanger -o {output.fq1o} -p {output.fq2o} -s {output.fqSo} -q 20 -l 20 -x -n 1> {log} ;"

# index rRNA sequence file for these genomes for read removal by bowtie2
rule refGenome_index: 
	input:
		fasta=RIBOSOME_RNA_fasta,
	params:
		"bowtie2idx/VE303_rRNA_bowtie2",
	output:
		bowtie2idx="bowtie2idx/VE303_rRNA_bowtie2.1.bt2",
	conda:
		"envs/bowtie2.yaml"
	shell:
		"bowtie2-build -q {input.fasta} {params} "
        
# remove reads that map to rRNA (computational depletion) 
# apparently this can be done just fine without trimming reads first...
rule bowtie2:
    input: 
        fq1="tmp/{sampleID}_filt/filt1.fq.gz", 
        fq2="tmp/{sampleID}_filt/filt2.fq.gz",
        bowtie2idx=rules.refGenome_index.output.bowtie2idx # put here, so rule bowtie2 only executed after rule refGenome_index done
    params:
        refGenome="bowtie2idx/VE303_rRNA_bowtie2",
        temp_prefix="1-bowtie2/{sampleID}_ref_ribosomes_unaligned.fastq", # just a prefix. 
    output:
        fqU1="1-bowtie2/{sampleID}_ref_ribosomes_unaligned.1.fastq",
        fqU2="1-bowtie2/{sampleID}_ref_ribosomes_unaligned.2.fastq",
    log:
        log="logs/bowtie2_{sampleID}_ref_ribosomes.txt",
    conda:
        "envs/bowtie2.yaml"
    shell:
        # 8 threads coded into json
        "bowtie2 --threads 8 -X 2000 --no-mixed --dovetail -x {params.refGenome} -1 {input.fq1} -2 {input.fq2} --un-conc-gz {params.temp_prefix} 2> {log} "

# run bowtie qc to get an idea of how well wet-lab rRNA depletion worked
rule bowtie2qc:
    input:
        bowtie2_logs = expand("logs/bowtie2_{sampleID}_ref_ribosomes.txt", sampleID=SAMPLE_ls),
    params:
        ref_str="ribosomes",
        out_file_root="1-bowtie2/alignment_stats",
    output:
        alignment_stats = "1-bowtie2/alignment_stats.csv",
    conda:
        "bowtie2QC",
    shell:
        "python3 {PROJECT_SCRIPTS_DIRECTORY}/bowtie2qc.py -s {spls} -d {CURRENT_DIRECTORY} -r {params.ref_str} -o {params.out_file_root}"


#build transcriptome index file using kallisto
rule kallistoidx:
    input:
        genome=ANNOT_REF_GENOME,
    output:
        index="VE303-04_kallisto.idx",
    conda:
        "kallisto",
    shell:
        "kallisto index -i {output.index} {input.genome}"

#use unaligned reads as input for kallisto
rule kallistoq_paired:
    input:
        fqU1="1-bowtie2/{sampleID}_ref_ribosomes_unaligned.1.fastq",
        fqU2="1-bowtie2/{sampleID}_ref_ribosomes_unaligned.2.fastq",
        index=rules.kallistoidx.output.index,
    output:
        "2-kallisto/{sampleID}/abundance.h5",
        "2-kallisto/{sampleID}/abundance.tsv",
        "2-kallisto/{sampleID}/run_info.json",
    params:
        threads=1
    conda:
        "kallisto",
    shell:
        "kallisto quant -i {input.index} -o 2-kallisto/{wildcards.sampleID} -b 100 -l 75 -s 2 {input.fqU1} {input.fqU2}"
        
# This rule assumes you're collecting all `abundance.tsv` outputs and running tximport from R
    # note: assuming 1 transcript per gene (appropriate for bacteria)
    # from ChatGPT -- need to validate
rule tximport_deseq2:
    input:
        abundance=expand("2-kallisto/{sampleID}/abundance.tsv", sampleID=SAMPLE_ls),
        metadata="metadata.csv"
    output:
        dds_rds="3-DESeq2/dds.rds",
        results_csv="3-DESeq2/deseq_results.csv",
        interp_csv="3-DESeq2/deseq2_interpretation.csv"
    conda:
        "deseq2",
    script:
        "scripts/tximport_deseq2.R"

        
#pseudoalign reads to transcriptome using kallisto
# rule kallistoq_single:
# 	input:
# 		filtfq1="/home/markey/mit_lieberman/projects/laura/2023_cacneso2_rnaseq/3-bowtie2/{sampleID}_ref_ribosomes_unaligned.1.fastq",
# 		index=rules.kallistoidx.output.index,
# 	output:
# 		"k_out/{sampleID}_filt_ribosomes/abundance.h5",
# 		"k_out/{sampleID}_filt_ribosomes/abundance.tsv",
# 		"k_out/{sampleID}_filt_ribosomes/run_info.json",
# 	params:
# 		threads=1
# 	shell:
# 		"kallisto quant -i {input.index} -o k_out/{wildcards.sampleID}_filt_ribosomes -b 100 --single -l 75 -s 2 {input.filtfq1}"



