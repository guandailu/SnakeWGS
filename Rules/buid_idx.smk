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
