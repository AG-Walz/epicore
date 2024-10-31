import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from matplotlib.ticker import MaxNLocator

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


def get_consensus_epitopes(protein_df, min_epi_len=9):
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


def gen_epitope(input_tsv, proteome_file, min_overlap=11, max_step_size=5,min_epi_len=9):
    '''
     input:
        - input_tsv: MHCquant output
        - protein_df: pandas DataFrame, one protein per row
        - min_overlap: minimum overlap length of two epitopes for a consensus epitope
        - max_step_size: maximum distance between the start positions of two epitopes for a consensus epitope
    output:
        - for each protein: list of core and whole peptides
    '''
    # get all proteins and their associated peptides
    protein_df = prot_pep_link(input_tsv)

    protein_df['end'] = protein_df.apply(lambda row: [x for _, x in sorted(zip(row['start'], row['end']), key=lambda y: y[0])],axis=1) 
    protein_df['sequence'] = protein_df.apply(lambda row: [x for _, x in sorted(zip(row['start'], row['sequence']), key=lambda y: y[0])],axis=1) 
    protein_df['start'] = protein_df['start'].apply(lambda x: sorted([int(i) for i in x])) 

    # group peptides 
    protein_df = group_peptides(protein_df, min_overlap, max_step_size, min_epi_len)
    protein_df = gen_landscape(protein_df,proteome_file)
    protein_df = get_consensus_epitopes(protein_df, min_epi_len=9)
    return protein_df


def map_pep_core(input_tsv, protein_df):
    '''
    input:
        - input_tsv: MHCquant output
        - protein_df: pandas DataFrame with one protein per row and all core and whole epitopes matched to that position 
    output:
        - input_tsv with two more columns, one for the whole epitope and one for the core epitope
    '''
    
    MHCquant_out = pd.read_csv(input_tsv, delimiter='\t')
    MHCquant_out['whole_epitopes'] = [[] for _ in range(len(MHCquant_out))]
    MHCquant_out['core_epitopes'] = [[] for _ in range(len(MHCquant_out))]
    
    for r, row in MHCquant_out.iterrows():

        # get protein 
        for mapping, accession in enumerate(row['accessions'].split(';')):

            # get start and end index for that group 
            start = row['start'].split(';')[mapping]
            end = row['end'].split(';')[mapping]
            sequence = row['sequence']

            prot_row = protein_df[(protein_df['accession'] == accession) & protein_df['sequence'].map(lambda x: sequence in x)]
            idx = [i for i, x in enumerate(zip(prot_row['start'].to_list()[0],prot_row['end'].to_list()[0])) if (x[0] == int(start) and x[1] == int(end))]
            if len(idx) > 1:
                #check if multiple occurrence due to modification
                wo_mod = [re.sub(r"\(.*?\)","",prot_row['sequence'].to_list()[0][i]) for i in idx]
                if len(set(wo_mod)) > 1:
                    print('Something went wrong with the mapping!')
                  

            mapped_group = prot_row['sequence_group_mapping'].to_list()[0][idx[0]]
            MHCquant_out.at[r,'core_epitopes'].append(prot_row['consensus_epitopes'].to_list()[0][mapped_group])
            MHCquant_out.at[r,'whole_epitopes'].append(prot_row['whole_epitopes'].to_list()[0][mapped_group])
    
        # convert list to ; separated strings
        MHCquant_out.at[r,'core_epitopes'] = ';'.join(MHCquant_out.at[r,'core_epitopes'])
        MHCquant_out.at[r,'whole_epitopes'] = ';'.join(MHCquant_out.at[r,'whole_epitopes'])
    return MHCquant_out


def vis_prot(protein_df, accession, proteome, plot_path=''):
    '''
    input:
        - protein_df: pandas DataFrame with one protein per row and all core and whole epitopes matched 
        - accession: accession of protein for visualization
        - plot_path: location where plot gets saved  
    '''
    prot_row = protein_df[(protein_df['accession'] == accession)]
    prot_seq = get_prot_seg(accession, proteome)

    prot_landscape = [0 for aa in prot_seq]

    fig, ax = plt.subplots(figsize=(6, 2), layout='constrained')
    ax.yaxis.set_major_locator(tck.MultipleLocator())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=20))
    group_landscapes = []
    group_cores = []
    max_y = 0
    for group, landscape in enumerate(prot_row['landscape'].to_list()[0]):
        group_landscape = [0 for aa in prot_seq]
        group_core = [0 for aa in prot_seq]

        #get start position
        group_start = min(prot_row['grouped_peptides_start'].to_list()[0][group])
        for idx, position in enumerate(landscape):
            group_landscape[group_start+int(idx)] += position
        for pos in range(prot_row['core_epitopes_start'].to_list()[0][group],prot_row['core_epitopes_end'].to_list()[0][group]):
            group_core[pos] += 1

        group_landscapes.append(group_landscape)
        group_cores.append(group_core)

        if max(group_landscape) > max_y:
            max_y = max(group_landscape)
    
    rgba_red = [1,0,0,1]
    rgba_red_a = [1,0,0,0.5]
    rgba_blue = [0,0,1,1]
    rgba_blue_a = [0,0,1,0.5]
    for i in range(len(group_landscapes)):
        if i % 2 == 0:
            ax.bar(range(len(prot_landscape)),group_landscapes[i],width=1, color=[rgba_red if pos == 1 else rgba_red_a for pos in group_cores[i]])
        else:
            ax.bar(range(len(prot_landscape)),group_landscapes[i],width=1, color=[rgba_blue if pos == 1 else rgba_blue_a for pos in group_cores[i]])
    
    ax.set_title('Number of peptides mapped to each amino acid position and core epitopes of protein {}'.format(accession))
    ax.set_xlabel('Position in protein {}'.format(accession))
    ax.set_ylabel('Number of mapped peptides')
    print(max_y)
    plt.show()


def __main__():
    args = parse_input()
    mhcquant_out = args.input_tsv
    fasta_proteome = args.proteome
   
    protein_df = gen_epitope(mhcquant_out,fasta_proteome,11,5,9)
    out_linked = map_pep_core(mhcquant_out,protein_df)
    protein_df.to_csv('plateau_result.csv')
    out_linked.to_csv('mhcquant_link_groups.csv')
    vis_prot(protein_df,'sp|P10909|CLUS_HUMAN',fasta_proteome)


if __name__ == "__main__":
    __main__()

# TODO: check filtering of short peptides
