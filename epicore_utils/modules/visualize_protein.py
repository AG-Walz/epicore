"""
Visualizes the landscape of a protein.
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import MaxNLocator
import pandas as pd 
import numpy as np
import logging
logger = logging.getLogger(__name__)


def plot_peptide_length_dist(first_df: pd.DataFrame, second_df: pd.DataFrame, first_explode: str, second_explode: str, first_column: str, second_column: str, first_label: str, second_label: str) -> tuple[plt.figure,int,int]:
    """Visualize the distribution of lengths of sequences.
    
    Args:
        first_df: A pandas dataframe.
        second_df: A pandas dataframe.
        first_explode: The header of the column, that is used to explode the 
            dataframe. 
        second_column: The header of the column, that is used to explode the 
            dataframe. 
        first_column: The header of the column, that holds the sequences in 
            first_df, for which the length distribution will be plotted. 
        second_column: The header of the column, which holds the sequences in 
            second_df, for which the length distribution will be plotted.  
        first_label: Label for the values of first_column in the plot.
        second_label: Label for the values of second_column in the plot. 

    Returns:
        A tuple including a matplotlib figure and two integers. The figure is a 
        histogram visualizing the length distribution of peptides and epitopes. 
        The two integers are the number of peptides and the number of epitopes. 
    """
    
    first_long = first_df.explode(first_explode)
    second_long = second_df.explode(second_explode)
    
    logger.info(f'{len(first_long)} peptides were reduced to {len(second_long)} epitopes.')

    # compute a histogram of the sequence lengths in first_column and second_column
    fig, ax = plt.subplots(layout='constrained')
    seq_first = first_long[first_column]
    seq_second = second_long[second_column]

    first_len = seq_first.map(lambda pep: len(pep)).to_list()
    second_len = seq_second.map(lambda pep: len(pep)).to_list()
    
    ax.hist(first_len, bins=np.arange(5,50,1), color='grey', label=first_label, alpha=0.6)
    ax.hist(second_len, bins=np.arange(5,50,1), color='red', label=second_label, alpha=0.6)

    ax.legend()
    ax.set_xlabel('length')
    ax.set_ylabel('count')
    return fig, len(first_long), len(second_long)


def plot_protein_landscape(protein_df: pd.DataFrame, accession: str, proteome_dict: dict[str,str]) -> plt.figure:
    """Visualize the landscape of a protein.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        accession: The string of a protein accession.
        proteome_dict: A dictionary containing the reference proteome.

    Returns:
        A matplotlib bar plot that visualizes the peptide and core epitope 
        distribution across the sequence of the protein with the provided 
        accession. 
    """

    prot_row = protein_df[(protein_df['accession'] == accession)]
    if len(prot_row) == 0:
        raise Exception('The accession {} is not in your input data.'.format(accession))

    prot_seq = proteome_dict[accession]

    prot_landscape = [0 for _ in prot_seq]
    
    fig_width = max(1,round(len(prot_landscape)/50))
    fig_height = 3

    max_intens = 0

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), layout='constrained')

    for group, landscape in enumerate(prot_row['landscape'].iloc[0]):

        if group % 3 == 0:
            color = 'red'
        elif group % 3 == 1:
            color = 'blue'
        else: 
            color = 'green'

        group_start = min(prot_row['grouped_peptides_start'].iloc[0][group])
        for idx, position in enumerate(landscape):
            ax.bar(group_start+int(idx),position,width=1, alpha=0.4, color=color)
        for pos in range(prot_row['core_epitopes_start'].iloc[0][group],prot_row['core_epitopes_end'].iloc[0][group]):
            ax.bar(pos,0.5,width=1,color=color)

        max_intens = max(max_intens, max(landscape))

    ybins=max_intens/10
    ax.yaxis.set_major_locator(MaxNLocator(nbins=max(ybins, 5)))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=max(fig_width/10, 15)))

    ax.set_title('Number of peptides mapped to each amino acid position and core epitopes of protein {}'.format(accession))
    ax.set_xlabel('Position in protein {}'.format(accession))
    ax.set_ylabel('Number of mapped peptides')
    return fig



def plot_core_mapping_peptides_hist(epitope_df: pd.DataFrame) ->  plt.figure:
    """A histogram of the number of peptides mapped to each epitope. 
    
    Args:
        epitope_df: A dataframe containing one epitope and 
            information of that epitope per row. 

    Returns:
        A histogram visualizing the number of peptides mapped to 
        each epitope. 
    """
    fig, ax = plt.subplots(layout='constrained')
    n_peps = epitope_df['grouped_peptides_sequence'].apply(lambda sequences: len(sequences.split(',')))
    ax.hist(n_peps,bins=np.arange(1,max(n_peps)+1,1))
    ax.set_yscale('log')
    ax.set_xlabel('number of peptides contributing to epitope')
    ax.set_ylabel('count')
    return fig


