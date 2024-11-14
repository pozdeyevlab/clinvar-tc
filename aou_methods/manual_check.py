"""
In order to confirm that the published results are correct, the variants in the following genes were manually confirmed to be listed as 'pathogenic' or 'likely pathogenic' in clinvar. 

This was done by querying the variant level results file, filtering for each gene, and then searching each varaint in the clinvar data base. 

The code below loops through a list of significant genes, and prints the variants in list format.
"""
import pandas as pd
import numpy as np
import defopt
from typing import List

def main(*, variant_results_file: str, list_of_genes: List[str])-> None:
    """
    Helper for manually checking variants pathogenicity

    :param variant_results_file: Path to variant level results from regression script
    :param list_of_genes: List of genes to check
    """
    vars: pd.DataFrame = pd.read_csv(variant_results_file, sep = '\t', header = 0)
    for gene in list_of_genes:
        tmp = vars[vars.gene == gene]
        tmp.to_csv(f'{gene}_check.tsv', sep='\t')

if __name__ == '__main__':
    defopt.run(main)