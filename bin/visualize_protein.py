"""
Visualizes the landscape of a protein.
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import MaxNLocator
import pandas as pd 
import numpy as np
import webbrowser 
import ast
import logging
logger = logging.getLogger(__name__)


def vis_pepdist(first_df: pd.DataFrame, second_df: pd.DataFrame, first_explode: str, second_explode: str, first_column: str, second_column: str, first_label: str, second_label: str) -> tuple[plt.figure,int,int]:
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


def vis_prot(protein_df: pd.DataFrame, accession: str, proteome_dict: dict[str,str], plot_path='') -> plt.figure:
    """Visualize the landscape of a protein.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        accession: The string of a protein accession.
        proteome_dict: A dictionary containing the reference proteome.
        plot_path: The location where the visualization gets saved to.

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



def pep_core_hist(epitope_df: pd.DataFrame) ->  plt.figure:
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
    ax.hist(n_peps,bins=np.arange(0,max(n_peps),1))
    ax.set_yscale('log')
    ax.set_xlabel('number of peptides contributing to epitope')
    ax.set_ylabel('count')
    return fig


def gen_report(length_distribution: str, intensity_hist: str, epitope_df: pd.DataFrame, peps: int, epitopes: int, n_removed_peps: int, min_overlap: int, max_step_size: int, min_epi_length: int, evidence_file: str):
    """ Generates a report including the most important information of the results. 

    Args:
        length_distribution: The path to the file containing the length 
            distribution of peptides and epitopes. 
        intensity_hist: The path to the file containing a histogram of the 
            number of peptides contributing to each core. 
        epitope_df: A dataframe containing one epitope and information of that  
            epitope per row. 
        peps: Number of peptides in the evidence file. 
        epitopes: Number of epitopes generated.
        n_removed_peps: Number of peptides removed from the evidence file due 
            to the absence of their accessions in the proteome.
        min_overlap: A user defined parameter. 
        max_step_size: A user defined parameter. 
        min_epi_length: A user defined parameter. 
        evidence_file: Path to the evidence file.

    Returns:
        Returns a html report summarizing some information of the localplateau 
        result. The report for example includes a figure of the peptide/epitope 
        length distribution, a histogram of the number of peptides mapped to 
        each epitope and the ten epitopes with the highest number of mapped 
        peptides. 
    """

    # compute mean of epitope core sequence
    len_core_epitopes = epitope_df['consensus_epitopes'].apply(lambda sequence: len(sequence))
    mean_core_length = sum(len_core_epitopes) / len(len_core_epitopes)

    # compute mean of epitope whole sequence
    len_whole_epitopes = epitope_df['whole_epitopes'].apply(lambda sequence: len(sequence))
    mean_whole_length = sum(len_whole_epitopes) / len(len_whole_epitopes)

    epitope_df_sort = epitope_df.sort_values(by='landscape', key=lambda landscapes: landscapes.apply(lambda landscape: max(ast.literal_eval(landscape))), ascending=False).head(10).to_html()

 
    report_f = open('report.html', 'w')
    html_content = f''' <html>
                        <head>
                            <style> 
                                .row {{
                                    display: flex;
                                    justify-content: space-between;
                                }}
                                .column {{
                                    flex: 1;
                                    margin: 0 10px;
                                }}
                                h1{{text-align: center;}}
                            </style>
                            <title>Localplateau report</title>
                        </head>
                        <body>
                            <h1>Localplateau report</h1>
                            <div class="row">
                                <div class="column"><p> The histogram shows the number of peptides/epitopes of a certain length. {peps} peptides were reduced to {epitopes} epitopes.</p><iframe src="{length_distribution}", style="width:600px; height:500px;" frameborder="0"></iframe></div>
                                <div class="column"><p> The histogram visualizes how many peptides contribute to the different epitopes.<br></p><iframe src="{intensity_hist}", style="width:600px; height:500px;" frameborder="0"></iframe></div>
                            </div>
                            <p> The average length of an epitope core sequence is {mean_core_length}. <br> The average length of an epitope sequence is {mean_whole_length}.</p><br>
                            <p>{n_removed_peps} peptides were removed, since they do not appear in the provided evidence file.</p>
                            <p> The above results were achieved with the following parameters:<br> - min_epi_length:{min_epi_length}<br> - min_overlap:{min_overlap}<br> - max_step_size:{max_step_size}<br>The input evidence file: {evidence_file}</p>
                            <p> The ten epitopes with the highest number of mapped peptides:<br>{epitope_df_sort}</p>
                            </body>
                        </html>'''
    report_f.write(html_content)
    report_f.close()
    webbrowser.open_new_tab('report.html')
