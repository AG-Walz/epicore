import os
import pandas as pd 
import ast
import yaml
import click
import logging
import numpy as np 
import matplotlib.pyplot as plt, mpld3
import warnings
import re

from . import __version__
from epicore_utils.modules.compute_cores import compute_consensus_epitopes
from epicore_utils.modules.map_result import map_pep_core, gen_epitope_df
from epicore_utils.modules.visualize_protein import plot_protein_landscape, plot_peptide_length_dist, plot_core_mapping_peptides_hist
from epicore_utils.modules.parse_input import parse_input, proteome_to_dict
from epicore_utils.modules.generate_report import gen_report

import logging

logger = logging.getLogger(__name__)

class InputParameter(object):
    """This class contains parameters necessary for the epicore script.

    Attributes:
        min_epi_length (int, optional): An integer of the minimal length of a 
            consensus epitope.
        min_overlap (int, optional): An integer of the minimal overlap between
            two epitopes to be grouped to the same consensus epitope.
        max_step_size (int, optional): An integer of the maximal distance 
            between the start position of two epitopes to be grouped to the same consensus epitope.
        seq_column (str, optional): The string of the header of the column 
            containing peptide sequence information in the evidence file.
        protacc_column (str, optional): The string of the header of the column  
            containing protein accession information in the evidence file.
        intensity_column (str, optional): The string of the header of the 
            column containing intensity information in the evidence file.
        delimiter (str, optional): The delimiter that separates multiple 
            entries in one column in the evidence file.
        mod_patter (str, optional): A comma separated string with delimiters
            for peptide modifications
        out_dir (str): A string of the directory, were all output files will be 
            saved.
        prot_accession (str, optional): A comma separated string containing the 
            protein accession, for which the protein landscape will be visualized.
        start_column (str, optional): The string of the header of the column 
            containing the start positions of peptides in proteins.
        end_column (str, optional): The string of the header of the column 
            containing the end position of peptides in proteins.

    """
    def __init__(self,reference_proteome=None, min_epi_length=None, min_overlap=None, max_step_size=None, seq_column=None, protacc_column=None, intensity_column=None, delimiter=None, mod_pattern=None, out_dir=None, prot_accession=None, start_column=None, end_column=None, report=None, html=None):
        self.min_epi_length = min_epi_length
        self.min_overlap = min_overlap
        self.max_step_size = max_step_size
        self.seq_column = seq_column
        self.protacc_column = protacc_column
        self.intensity_column = intensity_column
        self.delimiter = delimiter
        self.mod_pattern = mod_pattern
        self.out_dir = out_dir
        self.prot_accession = prot_accession
        self.start_column = start_column
        self.end_column = end_column
        self.report = report
        self.html = html
        self.proteome_dict = proteome_to_dict(reference_proteome)
        self.reference_proteome = reference_proteome

@click.version_option(__version__, "--version", "-V")

@click.group()
@click.option('--reference_proteome',type=click.Path(exists=True), required=True)
@click.option('--out_dir', type=click.Path(), required=True)
@click.pass_context
def main(ctx, reference_proteome, out_dir):
    ctx.obj = InputParameter(reference_proteome=reference_proteome, out_dir=out_dir)
    
@click.option('--min_epi_length', type=click.INT, default=11)
@click.option('--min_overlap', type=click.INT, default=11)
@click.option('--max_step_size', type=click.INT, required=5)
@click.option('--seq_column', type=click.STRING, required=True)
@click.option('--protacc_column', type=click.STRING, required=True)
@click.option('--intensity_column', type=click.STRING)
@click.option('--delimiter', type=click.STRING, required=True)
@click.option('--mod_pattern', type=click.STRING)
@click.option('--prot_accession', type=click.STRING)
@click.option('--start_column', type=click.STRING)
@click.option('--end_column', type=click.STRING)
@click.option('--report', is_flag=True)
@click.option('--html', is_flag=True)
@click.command()
@click.option('--evidence_file',type=click.Path(exists=True), required=True)
@click.pass_context
def generate_epicore_csv(ctx,evidence_file, min_epi_length, min_overlap, max_step_size, seq_column, protacc_column, intensity_column, delimiter, mod_pattern, prot_accession, start_column, end_column, report, html):
    ctx.obj = InputParameter(ctx.obj.reference_proteome, min_epi_length, min_overlap, max_step_size, seq_column, protacc_column, intensity_column, delimiter, mod_pattern, ctx.obj.out_dir, prot_accession, start_column, end_column, report, html)
    if not os.path.exists(ctx.obj.out_dir):
        os.mkdir(ctx.obj.out_dir)
    logging.basicConfig(filename=f'{ctx.obj.out_dir}/epicore.log', level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
        
    # ----------------------
    #    Parse input file
    # ----------------------
    protein_df, n_removed_peps = parse_input(evidence_file, ctx.obj.seq_column, ctx.obj.protacc_column, ctx.obj.intensity_column, ctx.obj.start_column, ctx.obj.end_column, ctx.obj.delimiter, ctx.obj.proteome_dict, ctx.obj.mod_pattern)
    os.makedirs(ctx.obj.out_dir,exist_ok=True)

    # ----------------------
    # compute core epitopes, map peptides to cores
    # ----------------------
    protein_df = compute_consensus_epitopes(protein_df, ctx.obj.min_overlap, ctx.obj.max_step_size, ctx.obj.min_epi_length, ctx.obj.intensity_column, ctx.obj.mod_pattern)
    protein_df.to_csv(f'{ctx.obj.out_dir}/epicore_result.csv')
    pep_cores_mapping = map_pep_core(evidence_file,protein_df,ctx.obj.seq_column,ctx.obj.protacc_column,ctx.obj.start_column,ctx.obj.end_column,ctx.obj.intensity_column,ctx.obj.delimiter,ctx.obj.mod_pattern, ctx.obj.proteome_dict)
    pep_cores_mapping.to_csv(f'{ctx.obj.out_dir}/pep_cores_mapping.csv')


    # ----------------------
    # Reformat data and generate multiple plots
    # ----------------------
    # generate file with one epitope in each row
    epitope_df = gen_epitope_df(protein_df)
    epitope_df.to_csv(f'{ctx.obj.out_dir}/epitopes.csv')

    # compute length distribution of peptides and epitopes
    ext = os.path.splitext(evidence_file)[1]
    if ext == '.csv':
        evidence_df = pd.read_csv(evidence_file, delimiter=',')
    elif ext == '.tsv':
        evidence_df = pd.read_csv(evidence_file, delimiter='\t')
    elif ext == '.xlsx':
        evidence_df = pd.read_excel(evidence_file)
    evidence_df[ctx.obj.protacc_column] = evidence_df[ctx.obj.protacc_column].apply(lambda accessions: accessions.split(ctx.obj.delimiter))

    fig = plot_core_mapping_peptides_hist(epitope_df)
    if ctx.obj.html:
        fig.savefig(f'{ctx.obj.out_dir}/epitope_intensity_hist.svg')
        with open(f'{ctx.obj.out_dir}/epitope_intensity_hist.svg', 'r') as svg_file:
            svg_content = svg_file.read()
        svg_content = re.sub(r'<\?xml[^>]+\?>', '', svg_content)
        svg_content = re.sub(r'<!DOCTYPE[^>]+>', '', svg_content)
        html = f'<!DOCTYPE html> <html> <body>{svg_content}</body></html>'
        with open(f'{ctx.obj.out_dir}/epitope_intensity_hist.html','w') as f:
            f.write(html)
    else:
        fig.savefig(f'{ctx.obj.out_dir}/epitope_intensity_hist.svg')
    fig, peps, epitopes = plot_peptide_length_dist(evidence_df, epitope_df, ctx.obj.seq_column, 'whole_epitopes', ctx.obj.seq_column, 'whole_epitopes', 'peptides', 'whole epitopes')
    if ctx.obj.html:
        fig.savefig(f'{ctx.obj.out_dir}/length_distributions.svg')
        with open(f'{ctx.obj.out_dir}/length_distributions.svg', 'r') as svg_file:
            svg_content = svg_file.read()
        svg_content = re.sub(r'<\?xml[^>]+\?>', '', svg_content)
        svg_content = re.sub(r'<!DOCTYPE[^>]+>', '', svg_content)
        html = f'<!DOCTYPE html> <html> <body>{svg_content}</body></html>'
        with open(f'{ctx.obj.out_dir}/length_distributions.html','w') as f:
            f.write(html)
    else:
        fig.savefig(f'{ctx.obj.out_dir}/length_distributions.svg')
    
    # summarize some results
    if ctx.obj.report:
        gen_report(f'./{ctx.obj.out_dir}/length_distributions.svg', f'{ctx.obj.out_dir}/epitope_intensity_hist.svg', epitope_df, peps, epitopes, n_removed_peps, ctx, evidence_file,  f'{ctx.obj.out_dir}/epicore_result.csv')


@click.command()
@click.option('--epicore_csv',type=click.Path(exists=True), required=True)
@click.option('--protacc', type=click.STRING, required=True)
@click.pass_context
def plot_landscape(ctx,epicore_csv, protacc):
    if not protacc:
        raise Exception('No protein accession was provided. Please provide a protein accession')
    for accession in protacc.split(','):

        # read in precomputed protein coverage and epitope cores.
        protein_df = pd.read_csv(epicore_csv)

        protein_df['grouped_peptides_start'] = protein_df['grouped_peptides_start'].apply(ast.literal_eval)
        protein_df['core_epitopes_start'] = protein_df['core_epitopes_start'].apply(lambda cell: eval(cell, {"np": np}))
        protein_df['core_epitopes_end'] = protein_df['core_epitopes_end'].apply(lambda cell: eval(cell, {"np": np}))
        protein_df['landscape'] = protein_df['landscape'].apply(ast.literal_eval)

        if accession is not None:
            fig = plot_protein_landscape(protein_df,accession,ctx.obj.proteome_dict)
            fig.savefig(f'{ctx.obj.out_dir}/{accession}.pdf',bbox_inches='tight')
            

main.add_command(generate_epicore_csv)
main.add_command(plot_landscape)

if __name__ == '__main__':
    main()
