rule prepde_counts:
    input:
        quant_gtfs = expand(
            f"{config['stringtie_dir']}/{{sample}}/{{sample}}.stringtie.gtf",
            sample=config["samples"].keys()
        )

    output:
        gene_counts = f"{config['prepde_dir']}/gene_count_matrix.csv",
        transcript_counts = f"{config['prepde_dir']}/transcript_count_matrix.csv"

    params:
        list_file = f"{config['prepde_dir']}/prepde_list.txt",
        script = config["prepde_script"],
        outdir = config["prepde_dir"]

    shell:
        """
        mkdir -p {params.outdir}

        rm -f {params.list_file}
        touch {params.list_file}

        for gtf in {input.quant_gtfs}; do
            base=$(basename "$gtf" .stringtie.gtf)
            echo "$base $gtf" >> {params.list_file}
        done

        python {params.script} \
            -i {params.list_file} \
            -g {output.gene_counts} \
            -t {output.transcript_counts}
        """