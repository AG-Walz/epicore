import argparse
import os 
import shutil
import pandas as pd 
import ast

from bin.compute_cores import gen_epitope
from bin.map_result import map_pep_core
from bin.visualize_protein import vis_prot

    
def __main__():

    parser = argparse.ArgumentParser(description='Input File and Parameters')
    parser.add_argument('-input_tsv', type=str,
                        help='PATH to the MHCquant output in tsv format.')
    parser.add_argument('-proteome', type=str, 
                        help='PATH to the fasta file of the proteome used for identification.')
    parser.add_argument('-min_epi_length', type=int, required=False, default=9,
                        help='Minimum epitope length. If no core is found for that length the core next longest core gets chosen.')
    parser.add_argument('-max_step_size', type=int, required=False, default=5,
                        help='The maximum step size between two epitopes to still be grouped together.')
    parser.add_argument('-min_overlap', type=int, required=False, default=5,
                        help='The minimum overlap between two epitopes to still be grouped together.')
    parser.add_argument('-prot_accession', type=str, required=False,
                        help='Accessions of proteins for which peptide coverage and core epitopes should be visualized. The input should be a string with comma separated accessions.')
    parser.add_argument('-out', type=str, required=True,
                        help='Path to where output tsvs and pdfs should be stored.')
    parser.add_argument('-plateau_csv', type=str, required=False,
                        help='Path to localplateau tsv. If set all input arguments other than -prot_accession, -out and -plateau_csv will be ignored and the peptide coverage and core epitopes of the given proteins will get visualized based on the plateau_csv.')
    args = parser.parse_args()

    mhcquant_out = args.input_tsv
    fasta_proteome = args.proteome
    min_epi_length = args.min_epi_length
    min_overlap = args.min_overlap
    max_step_size = args.max_step_size
    prot_accession = args.prot_accession
    out_dir = args.out 
    plateau_csv = args.plateau_csv

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
        protein_df = gen_epitope(mhcquant_out, min_overlap, max_step_size, min_epi_length)
        out_linked = map_pep_core(mhcquant_out,protein_df)
        protein_df.to_csv(out_dir + '/plateau_result_two.csv')
        out_linked.to_csv(out_dir + '/mhcquant_link_groups_two.csv')

        # visualize result - examples
        # class one 
        #vis_prot(protein_df,'sp|P10909|CLUS_HUMAN',fasta_proteome)
        #vis_prot(protein_df,'sp|P01024|CO3_HUMAN',fasta_proteome)
        #vis_prot(protein_df,'sp|P04114|APOB_HUMAN',fasta_proteome,'sp|P04114|APOB_HUMAN.pdf')
        # class two 
        #vis_prot(protein_df,'sp|P02671|FIBA_HUMAN',fasta_proteome,'two_sp|P02671|FIBA_HUMAN.pdf') ### look at position 532 - 539 !!!
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
            vis_prot(protein_df,accession,fasta_proteome,out_dir + '/' + accession + '.pdf')




if __name__ == "__main__":
    __main__()