import pandas as pd
import os.path
import subprocess

#### input configuration ####
configfile: "/group/zhougrp2/dguan/ChickenSV/Pipeline/config.yaml"
genome = config["genome"]
gtf  =   config["gtf"]
chr = list(range(1, 39)) + ["Z", "W", "MT"]
sample_tab=pd.read_csv(config["sample_tab"], header=0, sep = "\t")
sample=sample_tab["SampleID"].drop_duplicates().to_list()

### Create temporary directory
subprocess.run(['mkdir', '-p', 'Temp'])

rule all:
    input:
        expand('08_bqsr_bam/{sample}_bqsr.bam', sample=sample),
        expand('05_aligned_reads/{sample}.unaligned.bam', sample=sample),
        #expand("11_vcfs/Chr{chr}.vcf.gz", chr = chr)

#### Indexing reference #####
rule build_idx:
    input:
        genome
    output:
        idx = expand(genome+"{suffix}", suffix=[".amb", ".ann"]),
        dict = genome + ".dict"
    conda:
        "Envs/bwa.yaml"
    threads: 2
    shell:
        """
        samtools faidx {input}
        bwa index {input}
        samtools dict {input} > {output.dict}
        """

### Trimming adaptors and low quality reads
ruleorder: trim_paired_reads > trim_reads
rule trim_reads:
    input:
        fq="03_raw_reads/{sample}.fq.gz"
    output:
        reads = '04_trimmed_reads/{sample}_trimmed.fq.gz',
        report = '04_trimmed_reads/{sample}.fq.gz_trimming_report.txt'
    conda:
        'Envs/trimgalore.yaml'
    threads: 12
    shell:
        'trim_galore --output_dir 04_trimmed_reads -q {config[mapq]} --clip_R1 5 --three_prime_clip_R1 5 --gzip --fastqc --basename {wildcards.sample} --trim-n --cores {threads} {input}'

rule trim_paired_reads:
    input:
        fq1="03_raw_reads/{sample}_R1.fq.gz",
        fq2="03_raw_reads/{sample}_R2.fq.gz"
    output:
        reads1 = '04_trimmed_reads/{sample}_val_1.fq.gz',
        reads2 = '04_trimmed_reads/{sample}_val_2.fq.gz',
        report1 = '04_trimmed_reads/{sample}_R1.fq.gz_trimming_report.txt',
        report2 = '04_trimmed_reads/{sample}_R2.fq.gz_trimming_report.txt'
    conda:
        'Envs/trimgalore.yaml'
    threads: 12
    shell:
        'trim_galore --output_dir 04_trimmed_reads -q {config[mapq]} --clip_R1 5 --three_prime_clip_R1 5 --clip_R2 5 --three_prime_clip_R2 5 --gzip --fastqc --basename {wildcards.sample} --trim-n --paired --cores {threads} {input}'

### Reads mappping
ruleorder: bwa_mem_paired > bwa_mem
rule bwa_mem:
    input:
        index = rules.build_idx.output,
        reads = rules.trim_reads.output.reads
    output:
        bam = temp('05_aligned_reads/{sample}.aligned.bam')
    threads: 12
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"
    conda:
        'Envs/bwa.yaml'
    shell:
        "bwa mem -M -t {threads} -R '{params.rg}' {config[genome]} {input.reads} | samtools view -bS - > {output}"

rule bwa_mem_paired:
    input:
        rules.build_idx.output,
        reads1 = rules.trim_paired_reads.output.reads1,
        reads2 = rules.trim_paired_reads.output.reads2
    output:
        bam = temp('05_aligned_reads/{sample}.aligned.bam')
    threads: 24
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA"
    conda:
        'Envs/bwa.yaml'
    shell:
        "bwa mem -M -t {threads} -R '{params.rg}' {config[genome]} {input.reads1} {input.reads2} | samtools view -bS - > {output}"

rule sort_bam:
    input:
        lambda wildcards: '05_aligned_reads/{}.aligned.bam'.format(wildcards.sample)
    output:
        bam=temp('06_sorted_bam/{sample}.sorted.bam'),
        bai=temp('06_sorted_bam/{sample}.sorted.bam.bai')
    #conda:
    #    'Envs/samtools.yaml'
    threads: 12
    shell:
        """
        samtools sort -@ 12 -O bam -o {output.bam} {input}
        samtools index -@ 12 {output.bam}
        """
### Extract unmapped reads
rule extract_unalignment:
    input:
        bam='06_sorted_bam/{sample}.sorted.bam',
        bai='06_sorted_bam/{sample}.sorted.bam.bai'
    output:
        unalign='05_aligned_reads/{sample}.unaligned.bam'
    threads: 6
    shell:
        """
        samtools view -b -f 4 -@ 6 {input.bam} > {output}
        """

### Remove PCR duplicates
rule remove_deplicates:
    input:
        rules.sort_bam.output.bam
    output:
        dedup=temp('07_dedup_bam/{sample}.dedup.bam'),
        bai=temp('07_dedup_bam/{sample}.dedup.bai'),
        metrics='07_dedup_bam/{sample}.metrics.txt'
    params:
        java_opts="-Xmx24G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=Temp" # optional
    threads: 12
    conda:
        'Envs/gatk4.yaml'
    shell:
        """
        gatk --java-options "{params.java_opts}" MarkDuplicates -I {input} -M {output.metrics} -O {output.dedup} --COMPRESSION_LEVEL 9 --CREATE_INDEX true --REMOVE_DUPLICATES true --TMP_DIR ./Temp
        """

### Base calibration
rule run_basecalibration:
    input:
        bam='07_dedup_bam/{sample}.dedup.bam'
    output:
        bqsr='08_bqsr_bam/{sample}_bqsr.table'
    params:
        java_opts="-Xmx24G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=Temp"  # optional
    threads: 12
    conda:
        'Envs/gatk4.yaml'
    shell:
        """
        gatk --java-options "{params.java_opts}" BaseRecalibrator --spark-runner LOCAL -I {input} --known-sites {config[dbsnp]} -O {output} -R {config[genome]} --tmp-dir ./Temp
        """

rule apply_bqsr:
    input:
        bam='07_dedup_bam/{sample}.dedup.bam',
        bqsr='08_bqsr_bam/{sample}_bqsr.table'
    output:
        bqsr='08_bqsr_bam/{sample}_bqsr.bam'
    params:
        java_opts="-Xmx24G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=Temp"  # optional
    threads: 12
    conda:
        'Envs/gatk4.yaml'
    shell:
        """
        gatk --java-options "{params.java_opts}" ApplyBQSR --spark-runner LOCAL -I {input.bam} -R {config[genome]} --bqsr-recal-file {input.bqsr} -O {output} --create-output-bam-index true --create-output-bam-md5 true --tmp-dir ./Temp
        """

rule call_gvcf:
    input:
        bam='08_bqsr_bam/{sample}_bqsr.bam'
    output:
        gvcf='09_gvcfs/{sample}_chr{chr}.g.vcf.gz'
    params:
        chr = "{chr}",
        java_opts="-Xmx24G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=Temp"  # optional
    threads: 12
    conda:
        'Envs/gatk4.yaml'
    shell:
        """
        gatk --java-options "{params.java_opts}" HaplotypeCaller -R {config[genome]} -I {input.bam} -O {output} -ERC GVCF -L {params.chr} --native-pair-hmm-threads {threads} --create-output-variant-index true --tmp-dir ./Temp -G Standard -G AS_Standard
        """
rule genomics_db_import:
    input:
        gvcfs=expand("09_gvcfs/{sample}_chr{chr}.g.vcf.gz", sample=sample, chr=chr),
        genome = genome
    output:
        db=directory("10_genomics_db/Chr{chr}"),
    params:
        intervals="{chr}",
        db="10_genomics_db/Chr{chr}",
        extra=lambda wildcards, input: ' --variant '.join(input.gvcfs),
        java_opts="-Xmx24G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=Temp"  # optional
    conda:
        "Envs/gatk4.yaml"
    threads: 12
    resources:
        mem_gb=24,
    shell:
        """
        gatk --java-options "{params.java_opts}" GenomicsDBImport -R {input.genome} {input.gvcfs} --intervals {params.intervals} --genomicsdb-workspace-path {params.db} -L {params.intervals} --reader-threads 12 --batch-size 50 --tmp-dir Temp
        """

rule genotype_gvcfs:
    input:
        gvcf="10_genomics_db/Chr{chr}",
        genome=genome
    output:
        vcf="11_vcfs/Chr{chr}.vcf.gz"
    params:
        db="gendb://10_genomics_db/Chr{chr}",
        java_opts="-Xmx24G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=Temp"
    conda:
        "Envs/gatk4.yaml"
    threads: 12
    resources:
        mem_gb=24
    shell:
        """
        gatk --java-options "{params.java_opts}" GenotypeGVCFs -R {input.genome} -V {params.db} -O {output.vcf} --tmp-dir Temp
        """
