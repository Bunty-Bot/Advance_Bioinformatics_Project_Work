rule samtools_index:
    input:
        bam_file = f"{config['samtools_dir']}/{{sample}}/{{sample}}.sorted.bam"

    params:
        outdir = config["samtools_dir"]    

    output:
        bam_index_file = f"{config['samtools_dir']}/{{sample}}/{{sample}}.sorted.bam.bai"

    shell:
        """

        mkdir -p {params.outdir}/{wildcards.sample}

        samtools index -o {output.bam_index_file} {input.bam_file}

        """