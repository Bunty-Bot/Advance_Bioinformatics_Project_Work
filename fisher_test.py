#!/usr/bin/env python3
"""
Author: Chris Ambagtsheer (student number: 1216414)
Description: perform fisher's test over GO-term (like) table
Usage: python3 fisher_test.py pathway_table.csv \
        deseq2_data.tsv output blast_file

Arguments:
    fisher_test.py: this script
    pathway_table.tsv: tab separated table that contains;
         pathway IDs ,
         pathway descriptions
         gene IDs
    deseq2_table.tsv: tab separated deseq2 output table
    output: output file
    blast_file: file that contains the conversion identifiers to go from
        reference genome identifiers to plantcyc identifiers
"""

#import statements
from sys import argv
from pathlib import Path
from scipy.stats import fisher_exact, false_discovery_control

def parse_data(pathway_file:Path, first_col:str, second_col:str) -> dict:
    """Summary: read data and parse into dictionary

    Parameters:
        pathway_file (Path): file reference to input .csv table
        first_col (str): column to retrieve, becomes keys
        second_col (str): second column to retrieve, becomes values

    Return
        first columnal values as keys and second columnal values as values
         (dict)
    """

    if pathway_file.suffix == '.csv':
        sep = ','
    elif pathway_file.suffix in ['.tsv','.txt']:
        sep = '\t'
    else:
        raise 'Invalid file type, choose .csv, .tsv or .txt'

    pw_dict = {}
    first_line = True
    with open(pathway_file) as file:
        for line in file:
            sp_line = line.strip().split(sep)

            if not first_line and len(sp_line)-1 >= secondarg_idx:
                if sp_line[secondarg_idx]:
                    if sp_line[firstarg_idx] not in pw_dict:
                        pw_dict[sp_line[firstarg_idx]] = [] #does key exist?...
                    pw_dict[sp_line[firstarg_idx]] += [sp_line[secondarg_idx]]

            #identify column indices
            elif first_line:
                first_line = False
                firstarg_idx = sp_line.index(first_col)
                secondarg_idx = sp_line.index(second_col)

    return pw_dict

def deseq2_dict_converter(deseq2_dict:dict, blast_f:Path) -> dict:
    """Summary: convert pepperbase identifiers into plantcyc identifiers

    Parameters:
        deseq2_dict (dict): contains pepperbase identifiers as keys and
            p-value as values
        blast_f (Path): contains BLAST entry plantcyc substitutes for
            pepperbase identifiers

    Return:
        updated deseq2 dictionary, now with plantcyc identifiers
    """

    #make blast file dictionary with plantcyc IDs as values and pepperbase
    #   IDs as keys
    b_dict = {}
    with open(blast_f) as b_file:
        for line in b_file:
            sp_line = line.strip().split('\t')
            bgene_id = sp_line[0][:-5]
            plantcyc_id = sp_line[1]
            if bgene_id not in b_dict:
                b_dict[bgene_id] = plantcyc_id

    #Replace pepperbase IDs with plantcyc IDs
    new_ds2_dict = {}
    for dgene_id in deseq2_dict:
        if dgene_id in b_dict:
            new_ds2_dict[b_dict[dgene_id].upper()] = (
                float(deseq2_dict[dgene_id][0]))

    return new_ds2_dict


def get_counts(pw_dict:dict, ds2_dict:dict) -> tuple:
    """Summary: get total pathway total gene counts in and out of pathway

    Parameters:
        pw_dict (dict): pathway description as keys and gene ids as values
        ds2_dict (dict): gene id as keys and p-value as value

    Return:
        total significant genes, total insignificant genes and
            per pathway significant and insignificant genes
    """

    sig_pw_genes = 0
    insig_pw_genes = 0
    tot_sig_genes = 0
    tot_insig_genes = 0
    fisher_dict = {}

    #get all significant and insignificant gene counts
    for gene in ds2_dict:
        if ds2_dict[gene] < 0.05:
            tot_sig_genes += 1
        else:
            tot_insig_genes += 1

    #get the significant and insignificant gene count per pathway
    for pw in pw_dict:
        for pw_gene in pw_dict[pw]:

            if pw_gene in ds2_dict:
                if float(ds2_dict[pw_gene]) < 0.05:
                    sig_pw_genes += 1
                else:
                    insig_pw_genes += 1
            else:
                insig_pw_genes += 1

        fisher_dict[pw] = (sig_pw_genes, insig_pw_genes)
        sig_pw_genes = 0
        insig_pw_genes = 0

    return fisher_dict, tot_sig_genes, tot_insig_genes

def fishers_test(fisher_dict:dict, tot_sig_genes:int, tot_insig_genes:int) \
        -> dict:
    """Summary: perform a Fisher's exact test on the pathways

    Parameters:
        fisher_dict (dict): contains gene IDs as keys and significant
            pathway genes, insignificant pathway genes, total significant
            genes and total insignificant genes as values, all (int).
        tot_sig_genes (int): number of total significant genes
        tot_insig_genes (int): number tot total insignificant genes

    Yield:
        GO-term-table-like information in list format
    """

    pwp_dict = {}
    for gene in fisher_dict:
        sig_pw_genes, insig_pw_genes = fisher_dict[gene]
        sig_non_pw_genes = tot_sig_genes - sig_pw_genes
        insig_non_pw_genes = tot_insig_genes - insig_pw_genes
        fisher_table = [[sig_pw_genes, insig_pw_genes],
                        [sig_non_pw_genes, insig_non_pw_genes]]

        pval = fisher_exact(fisher_table, alternative = 'two-sided')
        pw_size = sig_pw_genes + insig_pw_genes
        gene_ratio = sig_pw_genes / pw_size

        pwp_dict[gene] = [tot_sig_genes, tot_insig_genes, pw_size,
                        sig_pw_genes, gene_ratio, float(pval.pvalue)]
    return pwp_dict

def bonf_corr(pwp_dict:dict) -> dict:
    """Summary: perfomr BH and bonferroni correction.

    Parameters:
        pwp_dict (dict): pathway dictionary with p-values and pathway data

    Return:
        updated pwp_dict with p-adjusted values
    """
    pval_list = []
    sorted_pw = sorted(list(pwp_dict))
    for pw in sorted_pw:
        pval_list += [pwp_dict[pw][-1]]

    adjusted = false_discovery_control(pval_list, method='bh')

    for idx, pw in enumerate(pwp_dict):
        pwp_dict[pw] += [float(adjusted[idx])]

    return pwp_dict

def make_file(pwpadj_dict:dict, output):
    """Summary: write significant pathway IDs and their values to file.

    Parameters:
        pwpadj_dict (dict): contains pathway identifier as key and pathway
            data as values.
        output (Path): file reference to output file.
    """

    with open(output, 'w') as file:
        file.write('description\ttotal_genes\tsignificant_genes\tpathway_size'
                   '\tsignificant_pathway_genes\tgene_ratio\tp\tpadj\n')
        for pw in pwpadj_dict:
            description = pw
            total_genes = pwpadj_dict[pw][0] + pwpadj_dict[pw][1]
            significant_genes = pwpadj_dict[pw][0]
            pathway_size = pwpadj_dict[pw][2]
            significant_pathway_genes = pwpadj_dict[pw][3]
            gene_ratio = pwpadj_dict[pw][4]
            p = pwpadj_dict[pw][5]
            padj = pwpadj_dict[pw][6]

            file.write(f'{description}\t{total_genes}\t{significant_genes}\t'
                       f'{pathway_size}\t{significant_pathway_genes}\t'
                       f'{gene_ratio}\t{p}\t{padj}\n')


def main():
    """Main function of the module."""

    #retrieve files from command line
    pathway_file = Path(argv[1])
    deseq2_results = Path(argv[2])
    output = Path(argv[3])
    blast_f = Path(argv[4])

    #make pathway IDs keys and gene IDs values
    pw_dict = parse_data(pathway_file,
                         'pathway_name',
                         'gene_id')

    ds2_dict = parse_data(deseq2_results,
                         'gene_id',
                         'padj')

    # convert identifiers to match with plantcyc
    conv_ds2_dict = deseq2_dict_converter(ds2_dict, blast_f)

    #get the count data needed for fisher tests
    fisher_dict, tot_sig_genes, tot_insig_genes = get_counts(pw_dict,
                                                             conv_ds2_dict)

    #perform fishers exact test
    pwp_dict = fishers_test(fisher_dict, tot_sig_genes, tot_insig_genes)

    #correct for multiple testing
    pwpadj_dict = bonf_corr(pwp_dict)

    #write data to file
    make_file(pwpadj_dict, output)

if __name__ == '__main__':
    main()