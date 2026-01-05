rule fastp:
    input:
        readF = lambda wc: config["samples"][wc.sample]["read1"],
        readR = lambda wc: config["samples"][wc.sample]["read2"]

    output:
        readF = f"{config['fastp_dir']}/{{sample}}/{{sample}}_1.trim.fastq",
        readR = f"{config['fastp_dir']}/{{sample}}/{{sample}}_2.trim.fastq",
        json  = f"{config['fastp_dir']}/{{sample}}/{{sample}}.json",
        html  = f"{config['fastp_dir']}/{{sample}}/{{sample}}.html"

    params:
        outdir = config["fastp_dir"]

    shell:
        """
        mkdir -p {params.outdir}/{wildcards.sample}

        fastp -i {input.readF} -I {input.readR} \
              -o {output.readF} -O {output.readR} \
              --json {output.json} --html {output.html} \
              --thread 4
        """