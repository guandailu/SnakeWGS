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
        '../Envs/gatk4.yaml'
    shell:
        """
        gatk --java-options "{params.java_opts}" BaseRecalibrator --spark-runner LOCAL -I {input} --known-sites {config[dbsnp]} -O {output} -R {config[genome]} --tmp-dir ./Temp
        """

rule apply_bqsr:
    input:
        bam='07_dedup_bam/{sample}.dedup.bam',
        bqsr='08_bqsr_bam/{sample}_bqsr.table'
    output:
        bqsr=temp('08_bqsr_bam/{sample}_bqsr.bam')
    params:
        java_opts="-Xmx24G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=Temp"  # optional
    threads: 12
    conda:
        '../Envs/gatk4.yaml'
    shell:
        """
        gatk --java-options "{params.java_opts}" ApplyBQSR --spark-runner LOCAL -I {input.bam} -R {config[genome]} --bqsr-recal-file {input.bqsr} -O {output} --create-output-bam-index true --create-output-bam-md5 true --tmp-dir ./Temp
        """

rule call_gvcf:
    input:
        bam='08_bqsr_bam/{sample}_bqsr.bam'
    output:
        gvcf='09_gvcfs/{sample}.chr{chr}.g.vcf.gz'
    params:
        chr = "{chr}",
        java_opts="-Xmx4G -XX:ParallelGCThreads=2" 
    threads: 2
    conda:
        '../Envs/gatk4.yaml'
    shell:
        """
        gatk --java-options "{params.java_opts}" HaplotypeCaller -R {config[genome]} -I {input.bam} -O {output} -ERC GVCF -L {params.chr} --native-pair-hmm-threads {threads} --create-output-variant-index true --tmp-dir {config[tempdir]}
        """

rule genomics_db_import:
    input:
        gvcfs=lambda wildcards: ["09_gvcfs/" + s + ".chr" + wildcards.chr + ".g.vcf.gz" for s in sample],
        genome = genome
    output:
        db=directory("10_genomics_db/chr{chr}"),
    params:
        intervals="{chr}",
        db="10_genomics_db/chr{chr}",
        gvcfs=lambda wildcards, input: ' --variant '.join(input.gvcfs),
        java_opts="-Xmx4G -XX:ParallelGCThreads=2"  # optional
    conda:
        "../Envs/gatk4.yaml"
    threads: 2
    resources:
        mem_gb=4,
    shell:
        """
        gatk --java-options "{params.java_opts}" GenomicsDBImport -R {input.genome} {params.gvcfs} --intervals {params.intervals} --genomicsdb-workspace-path {params.db} -L {params.intervals} --reader-threads 12 --batch-size 50 --tmp-dir {config[tempdir]}
        """

rule genotype_gvcfs:
    input:
        gvcf=rules.genomics_db_import.output,
        genome=genome
    output:
        vcf="11_called_vcfs/chr{chr}.vcf.gz"
    params:
        db="gendb://10_genomics_db/chr{chr}",
        java_opts="-Xmx4G -XX:ParallelGCThreads=2"
    conda:
        "../Envs/gatk4.yaml"
    threads: 2
    resources:
        mem_gb=4
    shell:
        """
        gatk --java-options "{params.java_opts}" GenotypeGVCFs -R {input.genome} -V {params.db} -O {output.vcf} --tmp-dir {config[tempdir]}
        """

