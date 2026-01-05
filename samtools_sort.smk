rule samtools_sort:
    input:
        sam_file = f"{config['hisat2_dir']}/{{sample}}/{{sample}}.sam"

    params:
        outdir = config["samtools_dir"]

    output:
        bam_file = f"{config['samtools_dir']}/{{sample}}/{{sample}}.sorted.bam"

    shell:
        """
        mkdir -p {params.outdir}/{wildcards.sample}

        samtools sort {input.sam_file} -o {output.bam_file}
        """