rule deseq2:
    input:
        metadata_table = config["metadata"],
        count_table = f"{config['prepde_dir']}/gene_count_matrix.csv"
    output:
        csv_file = f"{config['deseq2_dir']}/{{comp}}.csv"
    params:
        script = config["deseq2_script"]
    threads: 4
    shell:
        """
        python3 {params.script} \
            {input.metadata_table} \
            {input.count_table} \
            {wildcards.comp} \
            {output.csv_file}
        """