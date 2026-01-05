rule hisat2:
    input:

        fr = f"{config['fastp_dir']}/{{sample}}/{{sample}}_1.trim.fastq",
        rr = f"{config['fastp_dir']}/{{sample}}/{{sample}}_2.trim.fastq"

    params:
        genome = config["index_prefix"],
        outdir = config["hisat2_dir"]

    output:
        sam_file = f"{config['hisat2_dir']}/{{sample}}/{{sample}}.sam",
        summary  = f"{config['hisat2_dir']}/{{sample}}/{{sample}}_report.txt"

    shell:
        """
        mkdir -p {params.outdir}/{wildcards.sample}

        hisat2 -1 {input.fr} -2 {input.rr} -x {params.genome} -S {output.sam_file} -p 4 --summary-file {output.summary} --dta
        """