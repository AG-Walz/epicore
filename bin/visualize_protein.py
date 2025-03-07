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


def vis_pepdist(first_df: pd.DataFrame, second_df: pd.DataFrame, first_column: str, second_column: str, first_label: str, second_label: str) -> plt.figure:
    """Visualize the distribution of lengths of sequences.
    
    Args:
        first_df: A pandas dataframe.
        second_df: A pandas dataframe.
        first_column: The header of the column, that holds the sequences in 
            first_df, for which the length distribution will be plotted. 
        second_column: The header of the column, which holds the sequences in 
            second_df, for which the length distribution will be plotted.  
        first_label: Label for the values of first_column in the plot.
        second_label: Label for the values of second_column in the plot.                
        delimiter: The delimiter that separates multiple entries in one column 
            in the evidence file.
        
    """
    
    first_long = first_df.explode(first_column)
    second_long = second_df.explode(second_column)
    
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
    plt.show()
    return fig


def vis_prot(protein_df: pd.DataFrame, accession: str, proteome_dict: dict[str,str], plot_path='') -> plt.figure:
    """Visualize the landscape of a protein.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        accession: The string of a protein accession.
        proteome_dict: A dictionary containing the reference proteome.
        plot_path: The location where the visualization gets saved to.
    """

    prot_row = protein_df[(protein_df['accession'] == accession)]
    if len(prot_row) == 0:
        raise Exception('The accession {} is not in your input data.'.format(accession))

    prot_seq = proteome_dict[accession]

    prot_landscape = [0 for _ in prot_seq]
    
    fig_width = max(1,round(len(prot_landscape)/50))
    fig_height = 3


    fig, ax = plt.subplots(figsize=(fig_width, fig_height), layout='constrained')
    ax.yaxis.set_major_locator(tck.MultipleLocator())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=20))

    for group, landscape in enumerate(prot_row['landscape'].to_list()[0]):

        if group % 3 == 0:
            color = 'red'
        elif group % 3 == 1:
            color = 'blue'
        else: 
            color = 'green'

        group_start = min(prot_row['grouped_peptides_start'].to_list()[0][group])
        for idx, position in enumerate(landscape):
            ax.bar(group_start+int(idx),position,width=1, alpha=0.4, color=color)
        for pos in range(prot_row['core_epitopes_start'].to_list()[0][group],prot_row['core_epitopes_end'].to_list()[0][group]):
            ax.bar(pos,0.5,width=1,color=color)

    ax.set_title('Number of peptides mapped to each amino acid position and core epitopes of protein {}'.format(accession))
    ax.set_xlabel('Position in protein {}'.format(accession))
    ax.set_ylabel('Number of mapped peptides')
    plt.show()
    return fig

