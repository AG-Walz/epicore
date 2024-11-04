import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import MaxNLocator

from bin.compute_cores import get_prot_seg


def vis_prot(protein_df, accession, proteome, plot_path=''):
    '''
    input:
        - protein_df: pandas DataFrame with one protein per row and all core and whole epitopes matched 
        - accession: accession of protein for visualization
        - proteome: fasta file used for identification
        - plot_path: location where plot gets saved  
    '''
    prot_row = protein_df[(protein_df['accession'] == accession)]
    prot_seq = get_prot_seg(accession, proteome)

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


