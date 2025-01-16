"""
Visualizes the landscape of a protein.
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import MaxNLocator

from bin.parse_input import get_prot_seq

def vis_prot(protein_df, accession, proteome_df, plot_path=''):
    """Visualize the landscape of a protein.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        accession: The string of a protein accession.
        proteome_df: TODO
        plot_path: The location where the visualization gets saved to.
    """
    prot_row = protein_df[(protein_df['accession'] == accession)]
    if len(prot_row) == 0:
        raise Exception('The accession {} is not in your input data.'.format(accession))

    prot_seq = get_prot_seq(accession, proteome_df)

    prot_landscape = [0 for _ in prot_seq]
    
    fig_width = max(1,round(len(prot_landscape)/50))
    fig_height = 3


    _, ax = plt.subplots(figsize=(fig_width, fig_height), layout='constrained')
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
    if plot_path != '':
        plt.savefig(plot_path)
    plt.show()


