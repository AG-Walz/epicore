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

#from bin.parse_input import get_prot_seq

def vis_pepdist(protein_df: pd.DataFrame, evidence_file: str, protacc_column: str, delimiter: str):
    """Visualize the distribution of the epitopes length.
    
    Args:
        protein_df: A pandas dataframe containing one protein per row.
        evidence_file: The string of the path to the evidence file.
        protacc_column: The string of the header of the column containing 
            protein accession information in the evidence file.
        delimiter: The delimiter that separates multiple entries in one column 
            in the evidence file.
        
    """
    
    protein_df_long_whole = protein_df.explode('whole_epitopes')
    protein_df_long_core = protein_df.explode('consensus_epitopes')

    evidence_df = pd.read_csv(evidence_file,sep='\t')
    evidence_df[protacc_column] = evidence_df[protacc_column].apply(lambda accessions: accessions.split(delimiter))
    evidence_df_long = evidence_df.explode(protacc_column)
    
    logger.info(f'{len(evidence_df_long)} peptides were reduced to {len(protein_df_long_whole)} epitopes.')

    # compute a histogram of the peptides/core epitopes/whole epitopes lengths
    fig, ax = plt.subplots(layout='constrained')
    peps_before = evidence_df_long['sequence']
    core_after = protein_df_long_core['consensus_epitopes']
    whole_after = protein_df_long_whole['whole_epitopes']

    peps_len = peps_before.map(lambda pep: len(pep)).to_list()
    core_len = core_after.map(lambda pep: len(pep)).to_list()
    whole_len = whole_after.map(lambda pep: len(pep)).to_list()
    
    ax.hist(peps_len, bins=np.arange(5,50,1), color='grey', label='peptides', alpha=0.4)
    ax.hist(core_len, bins=np.arange(5,50,1), color='red', label='core_epitopes', alpha=0.7)
    ax.hist(whole_len,bins=np.arange(5,50,1), color='blue', label='whole epitopes', alpha=0.7)

    ax.legend()
    ax.set_xlabel('length')
    ax.set_ylabel('count')
    plt.show()
    return fig


def vis_prot(protein_df: pd.DataFrame, accession: str, proteome_dict: dict[str,str], plot_path=''):
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

