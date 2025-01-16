import argparse
import os 
import shutil
import pandas as pd 
import ast
import yaml

from bin.compute_cores import gen_epitope
from bin.map_result import map_pep_core
from bin.visualize_protein import vis_prot
from bin.parse_input import parse_input
from bin.parse_input import proteome_to_df

def __main__():

    parser = argparse.ArgumentParser(description='Input File and Parameters')
    parser.add_argument('-input_tsv', type=str,
                        help='PATH to the evidence file.')
    parser.add_argument('-proteome', type=str, 
                        help='PATH to the fasta file of the proteome used for identification.')
    parser.add_argument('-params_file', type=str, required=False, default=9,
                        help='PATH to the params.yaml file.')
    args = parser.parse_args()

    evidence_file = args.input_tsv
    fasta_proteome = args.proteome
    params_file = args.params_file

    # read parameters defined in yaml file 
    with open(params_file,'r') as yaml_file:
        params = yaml.safe_load(yaml_file)

    print(params)

    min_epi_length = params['parameters']['min_epi_length']
    min_overlap = params['parameters']['min_overlap']
    max_step_size = params['parameters']['max_step_size']
    seq_column = params['parameters']['seq_column']
    protacc_column = params['parameters']['protacc_column']
    intensity_column = params['parameters']['intensity_column']
    delimiter = params['parameters']['delimiter']
    mod_pattern = params['parameters']['mod_pattern']
    plateau_csv = params['parameters']['plateau_csv']
    out_dir = params['parameters']['out_dir']
    prot_accession = params['parameters']['prot_accession']
    start_column = params['parameters']['start_column']
    end_column = params['parameters']['end_column']

    # check if all input paths exist
    if not os.path.isfile(evidence_file):
        raise Exception('The given evidence file does not exist.')
    if not os.path.isfile(fasta_proteome):
        raise Exception('The given proteome file does not exist.')
    
    proteome_df = proteome_to_df(fasta_proteome)

    if plateau_csv is None:
        # parse input and compute start and end positions of peptides in proteins if search engine output does not provide position
        protein_df = parse_input(evidence_file, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, proteome_df, mod_pattern)
        os.makedirs(out_dir,exist_ok=True)
            
        
        # compute core epitopes and map peptides to cores
        protein_df = gen_epitope(protein_df, min_overlap, max_step_size, min_epi_length, intensity_column, mod_pattern)
        protein_df.to_csv(out_dir + '/plateau_result.csv')
        out_linked = map_pep_core(evidence_file,protein_df,seq_column,protacc_column,start_column,end_column,intensity_column,delimiter,mod_pattern)
        out_linked.to_csv(out_dir + '/evidence_link_groups.csv')

        if prot_accession is not None:
            for accession in prot_accession.split(','):
                vis_prot(protein_df,accession,proteome_df,out_dir + '/' + accession + '.pdf')
    
    else:
        for accession in prot_accession.split(','):
            # read in precomputed peptide coverage and epitope cores.
            protein_df = pd.read_csv(plateau_csv)
            protein_df['grouped_peptides_start'] = protein_df['grouped_peptides_start'].apply(ast.literal_eval)
            protein_df['core_epitopes_start'] = protein_df['core_epitopes_start'].apply(ast.literal_eval)
            protein_df['core_epitopes_end'] = protein_df['core_epitopes_end'].apply(ast.literal_eval)
            protein_df['landscape'] = protein_df['landscape'].apply(ast.literal_eval)
            if prot_accession is not None:
                for accession in prot_accession.split(','):
                    vis_prot(protein_df,accession,proteome_df,out_dir + '/' + accession + '.pdf')       
    




if __name__ == "__main__":
    __main__()
