rule deepvariant_gvcf:
    input:
        bam="08_bqsr_bam/{sample}_bqsr.bam",
        ref=genome
    output:
        vcf="09_gvcfs/{sample}.chr{chr}.vcf.gz",
        gvcf="09_gvcfs/{sample}.chr{chr}.g.vcf.gz"
    params:
        model="WGS",
        intermed_dir="09_gvcfs",
        chr="{chr}"
    threads: 1
    resources:
        mem_mb= 2
    shell:
        """
        module load singularity
        singularity run --nv -B /usr/lib/locale/:/usr/lib/locale/ docker://google/deepvariant:1.4.0 /opt/deepvariant/bin/run_deepvariant --model_type={params.model} --ref={input.ref} --reads={input.bam} --regions {params.chr} --output_vcf={output.vcf} --output_gvcf={output.gvcf} --intermediate_results_dir {params.intermed_dir} --num_shards={threads}
        """

rule gvcf_merge:
    input:
        gvcfs=lambda wildcards: ["09_gvcfs/" + s + ".chr" + wildcards.chr + ".g.vcf.gz" for s in sample],
    output:
        bcf="10_bcf/chr{chr}.bcf"
    threads: 12
    resources:
        mem_mb= 12000
    shell:
        """
        config["glnexuscli"] --config DeepVariantWGS --more-PL --mem-gbytes {resources.mem_mb} --threads {threads} {input} > {output}
        """

rule bcf2vcf:
    input:
        bcf="10_bcf/chr{chr}.bcf",
    output:
        vcf="11_called_vcfs/chr{chr}.vcf.gz",
    threads: 12
    resources:
        mem_mb= 12000
    conda:
        '../Envs/bcftools.yaml'
    shell:
        """
        bcftools view {input} | bgzip -@ {threads} -c > {output}
        """

