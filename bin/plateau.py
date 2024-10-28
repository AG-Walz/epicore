import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np
import re

def parse_input():
    parser = argparse.ArgumentParser(description='Input File and Parameters')
    parser.add_argument('-input_tsv', type=str,
                        help='PATH to the MHCquant output in tsv format.')
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
        if protein.id == peptide:
            return protein.seq


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


def remove_short_peptides(protein_df, min_epi_len=9):
    '''
    input:
        - protein_df: pandas DataFrame, one protein per row
        - min_epi_len: minimum length of an epitope
    output:
        - protein_df: input pandas DataFrame, without the epitopes shorter than min_epi_len
    '''

    for r,row in protein_df.iterrows():
        new_start = []
        new_end = []
        for i in range(len(row['start'])):
            if row['end'][i] - row['start'][i] > min_epi_len:
                new_start.append(row['start'][i])
                new_end.append(row['end'][i])
        protein_df.at[r,'start'] = new_start
        protein_df.at[r,'end'] = new_end

    return protein_df


def group_peptides(protein_df, min_overlap=11, max_step_size=5, min_epi_len=9):
    '''
    input:
        - protein_df: pandas DataFrame, one protein per row
        - min_overlap: minimum overlap length of two epitopes for a consensus epitope
        - max_step_size: maximum distance between the start positions of two epitopes for a consensus epitope
        - min_epi_len: minimum length of an epitope
    output:
        - input DataFrame, with the start and end positions and peptide sequences grouped into consensus epitopes
    '''

    protein_df['grouped_peptides_start'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_end'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_sequence'] = [[] for _ in range(len(protein_df))]

    for r,row in protein_df.iterrows():
        
        start_pos =  row['start']
        end_pos = row['end']
        sequences = row['sequence']
        grouped_peptides_start = []
        grouped_peptides_end = []
        grouped_peptides_sequence = []
        
        for i in range(len(start_pos)-1):

            grouped_peptides_start.append(start_pos[i])
            grouped_peptides_end.append(end_pos[i])
            grouped_peptides_sequence.append(sequences[i])

            step_size = start_pos[i+1] - start_pos[i]
            pep_length = end_pos[i] - start_pos[i]

            # create new peptide group after each jump
            if (step_size >= max_step_size) and (pep_length <= step_size + min_overlap):
                protein_df.at[r,'grouped_peptides_start'].append(grouped_peptides_start)
                protein_df.at[r,'grouped_peptides_end'].append(grouped_peptides_end)
                protein_df.at[r,'grouped_peptides_sequence'].append(grouped_peptides_sequence)
                grouped_peptides_end = []
                grouped_peptides_start = []
                grouped_peptides_sequence = []

        # special case for last peptide match of protein
        if len(grouped_peptides_end) == 0:
            protein_df.at[r,'grouped_peptides_start'].append([start_pos[-1]])
            protein_df.at[r,'grouped_peptides_end'].append([end_pos[-1]])
            protein_df.at[r,'grouped_peptides_sequence'].append([sequences[-1]])
        else:
            grouped_peptides_start.append(start_pos[-1])
            grouped_peptides_end.append(end_pos[-1])
            grouped_peptides_sequence.append(sequences[-1])
            protein_df.at[r,'grouped_peptides_start'].append(grouped_peptides_start)
            protein_df.at[r,'grouped_peptides_end'].append(grouped_peptides_end)
            protein_df.at[r,'grouped_peptides_sequence'].append(grouped_peptides_sequence)

    return protein_df
            
def gen_landscape(protein_df,proteome_file, min_overlap=11, max_step_size=5, min_epi_len=9):
    '''
    input:
        - protein_df: pandas DataFrame, one protein per row
        - min_overlap: minimum overlap length of two epitopes for a consensus epitope
        - max_step_size: maximum distance between the start positions of two epitopes for a consensus epitope
    output:
        - protein DataFrame + epitope landscape for each epitope group + whole sequence of each group
    '''

    protein_df['landscape'] = [[] for _ in range(len(protein_df))]
    protein_df['whole_epitopes'] = [[] for _ in range(len(protein_df))]

    for r, row in protein_df.iterrows():
        for group, [pep_group_start, pep_group_end] in enumerate(zip(row['grouped_peptides_start'], row['grouped_peptides_end'])):
            start_idx = pep_group_start[0]
            group_landscape = np.zeros(max(pep_group_end)+1-min(pep_group_start)) 
            for pep_start, pep_end in zip(pep_group_start, pep_group_end):
                for pos in range(pep_start, pep_end+1):
                    # only for identification, for quantification add intensity here
                    group_landscape[pos-start_idx] += 1

            protein_df.at[r,'landscape'].append(group_landscape)

            # build whole group sequence without proteome file
            # TODO: add case for grouped peptides with different modifications?
            consensus_seq = ''
            consensus_pos = []
            for sequence, sequence_pos in zip(row['grouped_peptides_sequence'][group], row['grouped_peptides_start'][group]):
                sequence_pos = [i for i in range(sequence_pos, sequence_pos + len(sequence))]
                for aa, aa_pos in zip(sequence, sequence_pos):
                    if aa_pos not in consensus_pos:
                        consensus_seq += aa
                        consensus_pos.append(aa_pos)
                    else:
                        if consensus_seq[consensus_pos.index(aa_pos)] != aa:
                            ('We have a problem here!')
            protein_df.at[r,'whole_epitopes'].append(consensus_seq)
       
            # build whole group sequence with proteome file
            # prot_seq = get_prot_seg(protein_df.at[r,'accession'], proteome_file)
            # protein_df.at[r,'whole_epitopes'].append(prot_seq[min(pep_group_start):max(pep_group_end)])  

    return protein_df


def get_consensus_epitopes(protein_df, min_epi_len=9):
    '''
    input:
        - protein_df: pandas DataFrame, one protein per row
    output:
        - input DataFrame, with the consensus epitope of each epitope group
    '''
    protein_df['consensus_epitopes'] = [[] for _ in range(len(protein_df))]
    for r, row in protein_df.iterrows():
        for group,landscape in enumerate(row['landscape']):
            
            # build consensus epitopes
            intens = np.unique(landscape)
            intens[::-1].sort()
        
            # find intensity for which consensus epitope is at least min_epi_len long
            for intensity in intens:

                Z = landscape < intensity

                # get lengths of peptide sequences with intensities above the current threshold
                seqs_idx = np.where(np.diff(np.hstack(([False],~Z,[False]))))[0].reshape(-1,2)
                
                # get length of longest peptide subsequences with current intensity
                ce_start_pos = seqs_idx[np.diff(seqs_idx, axis=1).argmax(),0]
                current_pep_length = np.diff(seqs_idx, axis=1).max()
                
                # check if min_epi_length is fulfilled for that sequence
                if current_pep_length >= min_epi_len:

                    # get position of epitope in protein sequences
                    pep_in_prot_start = ce_start_pos
                    pep_in_prot_end = pep_in_prot_start + current_pep_length

                    # get consensus epitopes
                    whole_epitope_wo_mod =re.sub(r"\(.*?\)","",protein_df.at[r,'whole_epitopes'][group])
                    protein_df.at[r,'consensus_epitopes'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                    break

def gen_epitope(input_tsv, proteome_file, min_overlap=11, max_step_size=5,min_epi_len=9):
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

    protein_df['end'] = protein_df.apply(lambda row: [x for _, x in sorted(zip(row['start'], row['end']))],axis=1) 
    protein_df['sequence'] = protein_df.apply(lambda row: [x for _, x in sorted(zip(row['start'], row['sequence']))],axis=1) 
    protein_df['start'] = protein_df['start'].apply(lambda x: sorted([int(i) for i in x])) 

    # group peptides 
    protein_df = group_peptides(protein_df, min_overlap, max_step_size, min_epi_len)
    protein_df = gen_landscape(protein_df,proteome_file)
    protein_df = get_consensus_epitopes(protein_df, min_epi_len=9)
    #protein_df = gen_landscape(protein_df, min_overlap, max_step_size)
    #protein_df = remove_short_peptides(protein_df, min_epi_len)
    return protein_df



def __main__():
    args = parse_input()
    mhcquant_out = args.input_tsv
    fasta_proteome = args.proteome
   
    protein_df = gen_epitope(mhcquant_out,fasta_proteome,11,5,9)
    protein_df.to_csv('plateau_result.csv')

if __name__ == "__main__":
    __main__()

# TODO: check filtering of short peptides