### Trimming adaptors and low quality reads
ruleorder: trim_paired_reads > trim_reads
rule trim_reads:
    input:
        fq="03_raw_reads/{sample}.fq.gz"
    output:
        reads = '04_trimmed_reads/{sample}_trimmed.fq.gz',
        report = '04_trimmed_reads/{sample}.fq.gz_trimming_report.txt'
    conda:
        '../Envs/trimgalore.yaml'
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
        '../Envs/trimgalore.yaml'
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
        '../Envs/bwa.yaml'
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
        '../Envs/bwa.yaml'
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
        dedup='07_dedup_bam/{sample}.dedup.bam',
        bai='07_dedup_bam/{sample}.dedup.bai',
        metrics='07_dedup_bam/{sample}.metrics.txt'
    params:
        java_opts="-Xmx24G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=Temp" # optional
    threads: 12
    conda:
        '../Envs/gatk4.yaml'
    shell:
        """
        gatk --java-options "{params.java_opts}" MarkDuplicates -I {input} -M {output.metrics} -O {output.dedup} --COMPRESSION_LEVEL 9 --CREATE_INDEX true --REMOVE_DUPLICATES true --TMP_DIR ./Temp
        """
