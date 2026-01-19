#!/usr/bin/env python3
"""
Author: Chris Ambagtsheer (student number: 1216414)
Description: perform differential expression analysis on a set of gene counts
             table with expressional data.
Usage: python3 DeSeq2.py metadata_table.csv counts_file.csv
       control_library treated_library output_file
Note: This version keeps full library names including replicates.
"""

from sys import argv
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd
import numpy as np

def load_data() -> tuple:
    """Load command line arguments."""
    metadata_file = Path(argv[1])
    counts_file = Path(argv[2])
    control = argv[3]
    treated = argv[4]
    output = Path(argv[5])
    return metadata_file, counts_file, control, treated, output

def parse_csv(filepath: Path, control: str, treated: str):
    """Yield fastq filename and full library name."""
    start = True
    with open(filepath) as file:
        for line in file:
            l_spl = line.strip().split(',')  # CSV file
            if start:
                l_spl = [x.strip() for x in l_spl]
                if l_spl[0].startswith('\ufeff'):
                    l_spl[0] = l_spl[0].replace('\ufeff', '')
                fastq_idx = l_spl.index('Run')
                libr_name_idx = l_spl.index('Library Name')
                start = False
            else:
                fastq = l_spl[fastq_idx]
                libr_name = l_spl[libr_name_idx]  # full name
                yield fastq, libr_name

def get_index_and_conditions(metadata_file: Path, control: str, treated: str) -> tuple:
    """Retrieve lists of fastq names and library names."""
    index = []
    conditions = []
    for fastq, libr_name in parse_csv(metadata_file, control, treated):
        index.append(fastq)
        conditions.append(libr_name)
    return index, conditions

def get_counts(counts_file: Path, index: list) -> pd.DataFrame:
    """Retrieve counts data."""
    raw_counts = pd.read_csv(counts_file, sep=',', index_col=0).T
    return raw_counts.loc[index]

def get_metadata(index: list, condition: list, control: str, treated: str) -> pd.DataFrame:
    """Return metadata dataframe for DESeq2."""
    meta = pd.DataFrame({'condition': condition}, index=index)
    # Convert to categorical safely
    meta['condition'] = meta['condition'].astype('category')
    meta['condition'] = meta['condition'].cat.reorder_categories(meta['condition'].unique(), ordered=True)
    return meta

def execute_deseq2(counts_df: pd.DataFrame, metadata_df: pd.DataFrame, control: str, treated: str) -> DeseqStats:
    """Run DESeq2 pipeline."""
    counts_df = counts_df.astype(int)

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata_df,
        refit_cooks=True
    )
    dds.fit_size_factors()
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    if dds.refit_cooks:
        dds.refit()

    # Map control and treated to indices in categories
    categories = metadata_df['condition'].cat.categories
    contrast = np.array([
        categories.get_loc(control),
        categories.get_loc(treated)
    ])

    ds = DeseqStats(dds,
                    contrast=contrast,
                    alpha=0.05,
                    cooks_filter=True,
                    independent_filter=True)

    ds.run_wald_test()
    if ds.cooks_filter:
        ds._cooks_filtering()
    if ds.independent_filter:
        ds._independent_filtering()
    else:
        ds._p_value_adjustment()
    ds.summary()

    return ds

def main():
    metadata_file, counts_file, control, treated, output = load_data()
    index, conditions = get_index_and_conditions(metadata_file, control, treated)
    counts_df = get_counts(counts_file, index)
    metadata_df = get_metadata(index, conditions, control, treated)
    ds = execute_deseq2(counts_df, metadata_df, control, treated)
    ds.results_df.to_csv(output)

if __name__ == '__main__':
    main()
