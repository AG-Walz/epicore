import argparse
import os 
import shutil
import pandas as pd 
import ast

from bin.params import *
from bin.compute_cores import gen_epitope
from bin.map_result import map_pep_core
from bin.visualize_protein import vis_prot
from bin.parse_input import parse_input

def __main__():
    
    pep_position = ''

    # check if all input paths exist
    if not os.path.isfile(mhcquant_out):
        raise Exception('The given evidence file does not exist.')
    if not os.path.isfile(fasta_proteome):
        raise Exception('The given proteome file does not exist.')

    # parse input and compute start and end positions of peptides in proteins if search engine output does not provide position
    protein_df = parse_input(mhcquant_out, seq_column, protacc_column, delimiter, fasta_proteome, mod_delimiter, pep_position)

    if plateau_csv is None:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        else:
            user_check = input('The output directory already exists and the results will be OVERWRITTEN. Enter "yes" to continue or "no" to stop.')
            if user_check == 'yes':
                shutil.rmtree(out_dir)
                os.mkdir(out_dir)
            else:
                exit()
        
        # compute core epitopes and map peptides to cores
        protein_df = gen_epitope(protein_df, min_overlap, max_step_size, min_epi_length)
        protein_df.to_csv(out_dir + '/plateau_result.csv')
        out_linked = map_pep_core(mhcquant_out,protein_df,delimiter)
        out_linked.to_csv(out_dir + '/evidence_link_groups.csv')
        
        # visualize result - examples
        # class one 
        #vis_prot(protein_df,'sp|P10909|CLUS_HUMAN',fasta_proteome)
        #vis_prot(protein_df,'sp|P01024|CO3_HUMAN',fasta_proteome)
        #vis_prot(protein_df,'sp|P04114|APOB_HUMAN',fasta_proteome,'sp|P04114|APOB_HUMAN.pdf')
        # class two 
        # vis_prot(protein_df,'sp|P02671|FIBA_HUMAN',fasta_proteome,'two_sp|P02671|FIBA_HUMAN.pdf') ### look at position 532 - 539 !!!
        #vis_prot(protein_df,'sp|P04114|APOB_HUMAN',fasta_proteome,'sp|P04114|APOB_HUMAN.pdf')

        if prot_accession is not None:
            for accession in prot_accession.split(','):
                vis_prot(protein_df,accession,fasta_proteome,out_dir + '/' + accession + '.pdf')
    
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
                    vis_prot(protein_df,accession,fasta_proteome,out_dir + '/' + accession + '.pdf')       
    




if __name__ == "__main__":
    __main__()



#TODO: check if index for positions is correct for plotting etc.

