from Bio import SeqIO
import numpy as np
import pandas as pd
import re


def get_prot_seg(accession, proteome_file):
    '''
    input:
        - accession (str)
        - reference proteome (fasta)
    output:
        - protein sequence corresponding to accession 
    '''
    proteome = SeqIO.parse(open(proteome_file),'fasta')
    for protein in proteome:
        if protein.id == accession:
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


def remove_short_peptides(protein_df, min_epi_len):
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


def group_peptides(protein_df, min_overlap, max_step_size):
    '''
    input:
        - protein_df: pandas DataFrame, one protein per row
        - min_overlap: minimum overlap length of two epitopes for a consensus epitope
        - max_step_size: maximum distance between the start positions of two epitopes for a consensus epitope
    output:
        - input DataFrame, with the start and end positions and peptide sequences grouped into consensus epitopes
    '''

    protein_df['grouped_peptides_start'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_end'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_sequence'] = [[] for _ in range(len(protein_df))]
    protein_df['sequence_group_mapping'] = [[] for _ in range(len(protein_df))]

    for r,row in protein_df.iterrows():
        
        start_pos =  row['start']
        end_pos = row['end']
        sequences = row['sequence']
        grouped_peptides_start = []
        grouped_peptides_end = []
        grouped_peptides_sequence = []
        n_jumps = 0
        mapping = []

        for i in range(len(start_pos)-1):

            grouped_peptides_start.append(start_pos[i])
            grouped_peptides_end.append(end_pos[i])
            grouped_peptides_sequence.append(sequences[i])

            step_size = start_pos[i+1] - start_pos[i]
            pep_length = end_pos[i] - start_pos[i]
            mapping.append(n_jumps)
            # create new peptide group after each jump
            if (step_size >= max_step_size) and (pep_length <= step_size + min_overlap):
                protein_df.at[r,'grouped_peptides_start'].append(grouped_peptides_start)
                protein_df.at[r,'grouped_peptides_end'].append(grouped_peptides_end)
                protein_df.at[r,'grouped_peptides_sequence'].append(grouped_peptides_sequence)
                n_jumps += 1
                grouped_peptides_end = []
                grouped_peptides_start = []
                grouped_peptides_sequence = []

        # special case for last peptide match of protein
        if len(grouped_peptides_end) == 0:
            protein_df.at[r,'grouped_peptides_start'].append([start_pos[-1]])
            protein_df.at[r,'grouped_peptides_end'].append([end_pos[-1]])
            protein_df.at[r,'grouped_peptides_sequence'].append([sequences[-1]])
            mapping.append(n_jumps)
        else:
            grouped_peptides_start.append(start_pos[-1])
            grouped_peptides_end.append(end_pos[-1])
            grouped_peptides_sequence.append(sequences[-1])
            protein_df.at[r,'grouped_peptides_start'].append(grouped_peptides_start)
            protein_df.at[r,'grouped_peptides_end'].append(grouped_peptides_end)
            protein_df.at[r,'grouped_peptides_sequence'].append(grouped_peptides_sequence)
            mapping.append(n_jumps)
        protein_df.at[r,'sequence_group_mapping'] = mapping

    return protein_df


def gen_landscape(protein_df):
    '''
    input:
        - protein_df: pandas DataFrame, one protein per row
    output:
        - protein DataFrame + epitope landscape for each epitope group + whole sequence of each group
    '''

    protein_df['landscape'] = [[] for _ in range(len(protein_df))]
    protein_df['whole_epitopes'] = [[] for _ in range(len(protein_df))]

    for r, row in protein_df.iterrows():
        for group, [pep_group_start, pep_group_end] in enumerate(zip(row['grouped_peptides_start'], row['grouped_peptides_end'])):
            start_idx = pep_group_start[0]
            group_landscape = [0 for _ in range(max(pep_group_end)+1-min(pep_group_start))]#np.zeros(max(pep_group_end)+1-min(pep_group_start)) 
            for pep_start, pep_end in zip(pep_group_start, pep_group_end):
                for pos in range(pep_start, pep_end+1):
                    # only for identification, for quantification add intensity here
                    group_landscape[pos-start_idx] += 1

            protein_df.at[r,'landscape'].append(group_landscape)

            # build whole group sequences
            consensus_seq = ''
            consensus_pos = []
            for sequence, sequence_pos in zip(row['grouped_peptides_sequence'][group], row['grouped_peptides_start'][group]):
                sequence = re.sub(r"\(.*?\)","",sequence)
                sequence_pos = [i for i in range(sequence_pos, sequence_pos + len(sequence))]
                for aa, aa_pos in zip(sequence, sequence_pos):
                    if aa_pos not in consensus_pos:
                        consensus_seq += aa
                        consensus_pos.append(aa_pos)
                    else:
                        if consensus_seq[consensus_pos.index(aa_pos)] != aa:
                            ('Something with the mapping went wrong here!')
            protein_df.at[r,'whole_epitopes'].append(consensus_seq) 

    return protein_df


def get_consensus_epitopes(protein_df, min_epi_len):
    '''
    input:
        - protein_df: pandas DataFrame, one protein per row
    output:
        - input DataFrame, with the consensus epitope of each epitope group
    '''
    protein_df['consensus_epitopes'] = [[] for _ in range(len(protein_df))]
    protein_df['core_epitopes_start'] = [[] for _ in range(len(protein_df))]
    protein_df['core_epitopes_end'] = [[] for _ in range(len(protein_df))]

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
                    whole_epitope_wo_mod = protein_df.at[r,'whole_epitopes'][group]
                    protein_df.at[r,'consensus_epitopes'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                    protein_df.at[r,'core_epitopes_start'].append(pep_in_prot_start+min(row['grouped_peptides_start'][group]))
                    protein_df.at[r,'core_epitopes_end'].append(pep_in_prot_end+min(row['grouped_peptides_start'][group])) 
                    break
                
                # if no core with length > min_epi_length
                if intensity == intens[-1]:
                    pep_in_prot_start = ce_start_pos
                    pep_in_prot_end = pep_in_prot_start + current_pep_length
                    protein_df.at[r,'consensus_epitopes'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                    protein_df.at[r,'core_epitopes_start'].append(pep_in_prot_start+min(row['grouped_peptides_start'][group]))
                    protein_df.at[r,'core_epitopes_end'].append(pep_in_prot_end+min(row['grouped_peptides_start'][group]))

    return protein_df


def gen_epitope(protein_df, min_overlap, max_step_size, min_epi_len):
    '''
     input:
        - input_tsv: MHCquant output
        - protein_df: pandas DataFrame, one protein per row
        - min_overlap: minimum overlap length of two epitopes for a consensus epitope
        - max_step_size: maximum distance between the start positions of two epitopes for a consensus epitope
    output:
        - for each protein: list of core and whole peptides
    '''

    protein_df['end'] = protein_df.apply(lambda row: [x for _, x in sorted(zip(row['start'], row['end']), key=lambda y: y[0])],axis=1) 
    protein_df['sequence'] = protein_df.apply(lambda row: [x for _, x in sorted(zip(row['start'], row['sequence']), key=lambda y: y[0])],axis=1) 
    protein_df['start'] = protein_df['start'].apply(lambda x: sorted([int(i) for i in x])) 

    # group peptides 
    protein_df = group_peptides(protein_df, min_overlap, max_step_size)
    protein_df = gen_landscape(protein_df)
    protein_df = get_consensus_epitopes(protein_df, min_epi_len)

    return protein_df