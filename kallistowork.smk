rule kallisto_index:
    input:
        genome = config["kallisto"]["genome_cds"]
    output:
        index = config["kallisto"]["index"]
    shell:
        """
        kallisto index -i {output.index} {input.genome}
        """


rule kallisto_quant:
    input:
        read1 = lambda wc: config["samples"][wc.sample]["read1"],
        read2 = lambda wc: config["samples"][wc.sample]["read2"],
        index = config["kallisto"]["index"]
    output:
        abundance = (
            config["kallisto"]["outdir"] + "/{sample}/abundance.tsv"
        )
    params:
        outdir = lambda wc: config["kallisto"]["outdir"] + f"/{wc.sample}"
    threads: 4
    shell:
        """
        kallisto quant \
            -i {input.index} \
            -o {params.outdir} \
            -t {threads} \
            {input.read1} {input.read2}
        """


rule kallisto_filter:
    input:
        abundance = config["kallisto"]["outdir"] + "/{sample}/abundance.tsv"
    output:
        filtered = config["kallisto"]["outdir"] + "/filtered/{sample}_filtered.tsv"
    shell:
        """
        mkdir -p $(dirname {output.filtered})
        awk 'NR==1 || $5 > 0' {input.abundance} > {output.filtered}
        """


rule kallisto_change_name:
    input:
        abundance = config["kallisto"]["outdir"] + "/{sample}/abundance.tsv"
    output:
        flat = config["kallisto"]["outdir"] + "/flat/{sample}_abundance.tsv"
    shell:
        """
        mkdir -p $(dirname {output.flat})
        cp {input.abundance} {output.flat}
        """