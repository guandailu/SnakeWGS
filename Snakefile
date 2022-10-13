import pandas as pd
import os.path
import subprocess

#### input configuration ####
configfile: "/group/zhougrp2/dguan/ChickenSV/Pipeline/config.yaml"
genome = config["genome"]
gtf  = config["gtf"]
chr = list(range(1, 39)) + ["Z", "W", "MT"]
sample_tab=pd.read_csv(config["sample_tab"], header=0, sep = "\t")
sample=sample_tab["SampleID"].drop_duplicates().to_list()

### Create temporary directory
subprocess.run(['mkdir', '-p', 'Temp'])

rule all:
    input:
        #expand('08_bqsr_bam/{sample}_bqsr.bam', sample=sample),
        expand('07_dedup_bam/{sample}.dedup.bam', sample=sample),
        expand('05_aligned_reads/{sample}.unaligned.bam', sample=sample),
        expand("11_called_vcfs/chr{chr}.vcf.gz", chr = chr)

include: "Rules/buid_idx.smk"
include: "Rules/reads_preprocess.smk"
if config["caller"] == "gatk":
    include: "Rules/gatk_call.smk"
elif config["caller"] == "deepvariant":
        include: "Rules/deepvariant_call.smk"
