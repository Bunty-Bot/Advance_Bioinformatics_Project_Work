rule stringtie_quantify:
    input:
        bam = f"{config['samtools_dir']}/{{sample}}/{{sample}}.sorted.bam",
        gtf = config["gff3"]

    output:
        quant_gtf = f"{config['stringtie_dir']}/{{sample}}/{{sample}}.stringtie.gtf",
        abundance = f"{config['stringtie_dir']}/{{sample}}/{{sample}}.abundance.txt"

    params:
        threads = 4,
        outdir = config["stringtie_dir"]

    shell:
        """
        mkdir -p {params.outdir}/{wildcards.sample}

        stringtie {input.bam} -p {params.threads} -G {input.gtf} -o {output.quant_gtf} -A {output.abundance} -e
        """