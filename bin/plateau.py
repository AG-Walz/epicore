import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np


def parse_input():
    parser = argparse.ArgumentParser(description='Input File and Parameters')
    parser.add_argument('-input_tsv', type=str,
                        help='PATH to the mhcquant output in tsv format.')
    parser.add_argument('-proteome', type=str, 
                        help='PATH to the fasta file of the proteome used for identification.')
    args = parser.parse_args()
    return args



def get_prot_seg(peptide, proteome_file):
    '''
    input:
        - accession (str)
        - reference proteome (fasta)
    output:
        - protein sequence corresponding to accession 
    '''
    proteome = SeqIO.parse(open(proteome_file),'fasta')
    for protein in proteome:
        print(protein.id)
        print(protein.seq)


def prot_pep_link(input_tsv):
    '''
    input:
        - input_tsv
    output:
        - pandas DataFrame
            > row for each protein accession in input:
                protein accession 
                peptides matched 
                start positions of peptide in protein
                end positions of peptide in protein
    '''

    # read MHCquant output to pd.DataFrame
    peptides = pd.read_csv(input_tsv, delimiter='\t')
    peptides = peptides[['accessions','sequence', 'start', 'end']]
    proteins = pd.DataFrame(columns=['accession', 'sequence', 'start','end'])

    for _, peptide in peptides.iterrows():
        # get all proteins associated with the peptide
        prot_accessions = peptide['accessions'].split(';')
        for i, prot_accession in enumerate(prot_accessions):
            if prot_accession in proteins['accession'].values:
                proteins.loc[proteins['accession'] == prot_accession, 'sequence'] = proteins.loc[proteins['accession'] == prot_accession, 'sequence'].apply(lambda x: x + [peptide['sequence']])
                proteins.loc[proteins['accession'] == prot_accession, 'start'] = proteins.loc[proteins['accession'] == prot_accession, 'start'].apply(lambda x: x + [int(pos) for pos in [peptide['start'].split(';')[i]]])
                proteins.loc[proteins['accession'] == prot_accession, 'end'] = proteins.loc[proteins['accession'] == prot_accession, 'end'].apply(lambda x: x + [int(pos) for pos in [peptide['end'].split(';')[i]]])
            else:
                protein_entry = {'accession':prot_accession, 'sequence':[peptide['sequence']], 'start':[int(pos) for pos in [peptide['start'].split(';')[i]]], 'end':[int(pos) for pos in [peptide['end'].split(';')[i]]]}
                proteins.loc[len(proteins)] = protein_entry

    return proteins


def gen_landscape(protein_df, min_overlap=11, max_step_size=5):
    '''
    input:
        - protein_df: pandas DataFrame, one protein per row
        - min_overlap: minimum overlap length of two epitopes for a consensus epitope
        - max_step_size: maximum distance between the start positions of two epitopes for a consensus epitope
    output:
        - protein DataFrame + epitope landscape for each protein
    '''

    protein_df['landscape'] = np.nan
    protein_df['end'] = protein_df.apply(lambda row: [x for _, x in sorted(zip(row['start'], row['end']))],axis=1) 
    protein_df['start'] = protein_df['start'].apply(lambda x: sorted([int(i) for i in x])) 
    return protein_df



def gen_epitope(input_tsv, min_overlap=11, max_step_size=5):
    '''
     input:
        - input_tsv: MHCquant output
        - protein_df: pandas DataFrame, one protein per row
        - min_overlap: minimum overlap length of two epitopes for a consensus epitope
        - max_step_size: maximum distance between the start positions of two epitopes for a consensus epitope
    output:
        - consensus epitopes
    '''

    # get all proteins and their associated peptides
    protein_df = prot_pep_link(input_tsv)
    protein_df = gen_landscape(protein_df, min_overlap, max_step_size)
    return protein_df



def __main__():
    args = parse_input()
    mhcquant_out = args.input_tsv
    fasta_proteome = args.proteome
    protein_df = gen_epitope(mhcquant_out)

if __name__ == "__main__":
    __main__()