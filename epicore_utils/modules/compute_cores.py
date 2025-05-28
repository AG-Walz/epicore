"""
Computes core epitopes, by grouping overlapping peptides, together, building a landscape for each group and identifying plateaus with a defined minimal length in each landscape. 
"""

import numpy as np
import re
import pandas as pd


def group_peptides(protein_df: pd.DataFrame, min_overlap: int, max_step_size: int, intensity_column: str) -> pd.DataFrame:
    """Group the peptides to consensus epitopes.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        min_overlap: An integer of the minimal overlap between two epitopes     
            to be grouped to the same consensus epitope.
        max_step_size: An integer of the maximal distance between the start 
            position of two epitopes to be grouped to the same consensus 
            epitope.
        intensity_column: The header of the column containing the intensities
            of the peptides.

    Returns:
        The protein_df with the five to seven additional columns: 
        grouped_peptides_start, grouped_peptides_end, 
        grouped_peptides_sequence, grouped peptides_intensity, 
        core_epitopes_intensity, relative_core_intensity and 
        sequence_group_mapping. The first five column contain the start
        position, end position, sequence and intensity of the peptides grouped
        to consensus epitopes. The column relative_core_intensity contains the 
        relative core intensities of the consensus core epitopes. The column 
        sequence group mapping contains for each peptide to which consensus 
        epitope it is grouped.
    """
    # start, end, sequence and intensity of peptides of one group grouped together
    protein_df['grouped_peptides_start'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_end'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_sequence'] = [[] for _ in range(len(protein_df))]
    if intensity_column:
        protein_df['grouped_peptides_intensity'] = [[] for _ in range(len(protein_df))]
        # for each peptide group the total and relative intensity of the entire group
        protein_df['core_epitopes_intensity'] = [[] for _ in range(len(protein_df))]
        protein_df['relative_core_intensity'] = [[] for _ in range(len(protein_df))]
    # for each peptide the index of its group
    protein_df['sequence_group_mapping'] = [[] for _ in range(len(protein_df))]

    total_file_intensity = 0
    for r,row in protein_df.iterrows():
        
        start_pos =  row['start']
        end_pos = row['end']
        sequences = row['sequence']
        if intensity_column:
            intensity = row['intensity']

        grouped_peptides_start = []
        grouped_peptides_end = []
        grouped_peptides_sequence = []
        if intensity_column:
            grouped_peptides_intensity = []
            core_intensity = 0
        n_jumps = 0
        mapping = []

        for i in range(len(start_pos)-1):

            grouped_peptides_start.append(int(start_pos[i]))
            grouped_peptides_end.append(int(end_pos[i]))
            grouped_peptides_sequence.append(sequences[i])
            if intensity_column:
                grouped_peptides_intensity.append(intensity[i])

            step_size = int(start_pos[i+1]) - int(start_pos[i])
            pep_length = int(end_pos[i]) - int(start_pos[i])
            mapping.append(n_jumps)
            if intensity_column:
                core_intensity += float(intensity[i])
                total_file_intensity += float(intensity[i])

            # create new peptide group after each jump
            if (step_size >= max_step_size) and (pep_length <= step_size + min_overlap):
                protein_df.at[r,'grouped_peptides_start'].append(grouped_peptides_start)
                protein_df.at[r,'grouped_peptides_end'].append(grouped_peptides_end)
                protein_df.at[r,'grouped_peptides_sequence'].append(grouped_peptides_sequence)
                if intensity_column:
                    protein_df.at[r,'grouped_peptides_intensity'].append(grouped_peptides_intensity)
                    protein_df.at[r,'core_epitopes_intensity'].append(core_intensity)

                n_jumps += 1
                grouped_peptides_end = []
                grouped_peptides_start = []
                grouped_peptides_sequence = []
                grouped_peptides_intensity = []
                if intensity_column:
                    core_intensity = 0

        # special case for last peptide match of protein
        if len(grouped_peptides_end) == 0:
            protein_df.at[r,'grouped_peptides_start'].append([int(start_pos[-1])])
            protein_df.at[r,'grouped_peptides_end'].append([int(end_pos[-1])])
            protein_df.at[r,'grouped_peptides_sequence'].append([sequences[-1]])
            if intensity_column:
                protein_df.at[r,'grouped_peptides_intensity'].append([intensity[-1]])
                protein_df.at[r,'core_epitopes_intensity'].append(intensity[-1])
                total_file_intensity += float(intensity[-1])
            mapping.append(n_jumps)
        else:
            grouped_peptides_start.append(int(start_pos[-1]))
            grouped_peptides_end.append(int(end_pos[-1]))
            grouped_peptides_sequence.append(sequences[-1])
            if intensity_column:
                grouped_peptides_intensity.append(intensity[-1])
                core_intensity += float(intensity[-1])
                total_file_intensity += float(intensity[-1])
                protein_df.at[r,'grouped_peptides_intensity'].append(grouped_peptides_intensity)
                protein_df.at[r,'core_epitopes_intensity'].append(core_intensity)
            protein_df.at[r,'grouped_peptides_start'].append(grouped_peptides_start)
            protein_df.at[r,'grouped_peptides_end'].append(grouped_peptides_end)
            protein_df.at[r,'grouped_peptides_sequence'].append(grouped_peptides_sequence)
            mapping.append(n_jumps)

        protein_df.at[r,'sequence_group_mapping'] = mapping
    if intensity_column:
        protein_df['relative_core_intensity'] = protein_df['core_epitopes_intensity'].apply(lambda x: [float(ints)/total_file_intensity for ints in x])
    if intensity_column:
        protein_df['core_epitopes_intensity_all'] = protein_df.apply(lambda row: [row['core_epitopes_intensity'][i] for i in row['sequence_group_mapping']],axis=1)
        protein_df['relative_core_intensity_all'] = protein_df.apply(lambda row: [row['relative_core_intensity'][i] for i in row['sequence_group_mapping']],axis=1)
    
    return protein_df


def gen_landscape(protein_df: pd.DataFrame, mod_pattern: str) -> pd.DataFrame:
    """Compute the landscape of consensus epitope groups.
    
    Args:
        protein_df: A pandas dataframe containing one protein per row.
        mod_pattern A comma separated string with delimiters for peptide
            modifications

    Returns:
        The protein_df with two additional columns (landscape, whole_epitopes).
        The column landscape contains the landscapes of the consensus epitopes 
        group. The height of the landscape of a consensus epitope group at a 
        position is the number of peptides that contain that position. The 
        column whole_epitopes contains the entire sequence of the consensus 
        epitope group.

    Raises:
        Exception: If the mappings are contradictory.
    """

    protein_df['landscape'] = [[] for _ in range(len(protein_df))]
    protein_df['whole_epitopes'] = [[] for _ in range(len(protein_df))]
    protein_df['whole_epitopes_all'] = [[] for _ in range(len(protein_df))]

    for r, row in protein_df.iterrows():
        #seen_pep={}
        for group, [pep_group_start, pep_group_end] in enumerate(zip(row['grouped_peptides_start'], row['grouped_peptides_end'])):
            seen_pep={}
            start_idx = pep_group_start[0]
            group_landscape = [0 for _ in range(max(pep_group_end)+1-min(pep_group_start))]#np.zeros(max(pep_group_end)+1-min(pep_group_start)) 
            
            for pep_idx, pep_pos in enumerate(zip(pep_group_start, pep_group_end)):
                pep_start = pep_pos[0]
                pep_end = pep_pos[1]

                current_seq = row['grouped_peptides_sequence'][group][pep_idx]
                pattern = r'\(.*?\)'
                current_seq = re.sub(pattern,"",current_seq)
                pattern = r'\[.*?\]'
                current_seq = re.sub(pattern,"",current_seq)
                if mod_pattern:
                    pattern = re.escape(mod_pattern.split(',')[0]) + r'.*?' + re.escape(mod_pattern.split(',')[1])
                    current_seq = re.sub(pattern,"",current_seq)

                if str(pep_start) + '-' + str(pep_end) in seen_pep:
                    seen_seq = seen_pep[str(pep_start) + '-' + str(pep_end)]
                    check_seen = [(seq == current_seq) for seq in seen_seq]
                    # current_seq not in list for repetitive regions (dont increase landscape value for repetitive regions)
                    seen_pep[str(pep_start) + '-' + str(pep_end)].append(current_seq)
                else: 
                    # if sequence and position not seen before add the sequence to the landscape
                    seen_pep[str(pep_start) + '-' + str(pep_end)] = [current_seq]
                    for pos in range(pep_start, pep_end+1):
                        group_landscape[pos-start_idx] += 1

            protein_df.at[r,'landscape'].append(group_landscape)

            # build whole group sequences
            consensus_seq = ''
            consensus_pos = []
            for sequence, sequence_pos in zip(row['grouped_peptides_sequence'][group], row['grouped_peptides_start'][group]):
                sequence = re.sub(r"\(.*?\)","",sequence)
                sequence = re.sub(r'\[.*?\]',"",sequence)
                if mod_pattern:
                    pattern = re.escape(mod_pattern.split(',')[0]) + r'.*?' + re.escape(mod_pattern.split(',')[1])
                    sequence = re.sub(pattern,"",sequence)
                sequence_pos = [i for i in range(sequence_pos, sequence_pos + len(sequence))]
                for aa, aa_pos in zip(sequence, sequence_pos):
                    if aa_pos not in consensus_pos:
                        consensus_seq += aa
                        consensus_pos.append(aa_pos)
                    else:
                        if consensus_seq[consensus_pos.index(aa_pos)] != aa:
                            raise Exception('Something with the mapping went wrong! If you provided the start and end column please try again without providing the column header for these columns.')
            for _ in row['grouped_peptides_sequence'][group]:
                protein_df.at[r,'whole_epitopes_all'].append(consensus_seq)
            protein_df.at[r,'whole_epitopes'].append(consensus_seq)
        
    return protein_df


def get_consensus_epitopes(protein_df: pd.DataFrame, min_epi_len: int) -> pd.DataFrame:
    """Compute the consensus epitope sequence of each consensus epitope group.
    
    Args:
        protein_df: A pandas dataframe containing one protein per row.
        min_epi_length: An integer of the minimal length of a consensus epitope.

    Returns:
        The protein_df with one additional column, that contains the consensus epitope sequence of each consensus epitope group.
    """
    protein_df['consensus_epitopes'] = [[] for _ in range(len(protein_df))]
    protein_df['consensus_epitopes_all'] = [[] for _ in range(len(protein_df))]
    protein_df['core_epitopes_start'] = [[] for _ in range(len(protein_df))]
    protein_df['core_epitopes_end'] = [[] for _ in range(len(protein_df))]
    protein_df['core_epitopes_start_all'] = [[] for _ in range(len(protein_df))]
    protein_df['core_epitopes_end_all'] = [[] for _ in range(len(protein_df))]

    for r, row in protein_df.iterrows():
        for group,landscape in enumerate(row['landscape']):
            
            # build consensus epitopes
            total_counts = np.unique(landscape)
            total_counts[::-1].sort()
        
            # find total coverage for which consensus epitope is at least min_epi_len long
            for total_count in total_counts:

                Z = landscape < total_count

                # get lengths of peptide sequences with coverage above the current threshold
                seqs_idx = np.where(np.diff(np.hstack(([False],~Z,[False]))))[0].reshape(-1,2)
                
                # get length of longest peptide subsequences with current count
                ce_start_pos = seqs_idx[np.diff(seqs_idx, axis=1).argmax(),0]
                current_pep_length = np.diff(seqs_idx, axis=1).max()
                
                # check if min_epi_length is fulfilled for that sequence
                if current_pep_length >= min_epi_len:

                    # get position of epitope in protein sequences
                    pep_in_prot_start = ce_start_pos
                    pep_in_prot_end = pep_in_prot_start + current_pep_length

                    # get consensus epitopes
                    whole_epitope_wo_mod = protein_df.at[r,'whole_epitopes'][group]
                    for _ in row['grouped_peptides_sequence'][group]:
                        protein_df.at[r,'consensus_epitopes_all'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                        protein_df.at[r,'core_epitopes_start_all'].append(pep_in_prot_start+min(row['grouped_peptides_start'][group]))
                        protein_df.at[r,'core_epitopes_end_all'].append(pep_in_prot_end+min(row['grouped_peptides_start'][group]) - 1) 
                    protein_df.at[r,'consensus_epitopes'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                    protein_df.at[r,'core_epitopes_start'].append(pep_in_prot_start+min(row['grouped_peptides_start'][group]))
                    protein_df.at[r,'core_epitopes_end'].append(pep_in_prot_end+min(row['grouped_peptides_start'][group]) - 1)
                    break
                
                # if no core with length > min_epi_length
                if total_count == total_counts[-1]:
                    pep_in_prot_start = ce_start_pos
                    pep_in_prot_end = pep_in_prot_start + current_pep_length
                    whole_epitope_wo_mod = protein_df.at[r,'whole_epitopes'][group]
                    for _ in row['grouped_peptides_sequence'][group]:
                        protein_df.at[r,'consensus_epitopes_all'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                        protein_df.at[r,'core_epitopes_start_all'].append(pep_in_prot_start+min(row['grouped_peptides_start'][group]))
                        protein_df.at[r,'core_epitopes_end_all'].append(pep_in_prot_end+min(row['grouped_peptides_start'][group]) - 1)
                    protein_df.at[r,'consensus_epitopes'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                    protein_df.at[r,'core_epitopes_start'].append(pep_in_prot_start+min(row['grouped_peptides_start'][group]))
                    protein_df.at[r,'core_epitopes_end'].append(pep_in_prot_end+min(row['grouped_peptides_start'][group]) - 1)
                    
    protein_df['proteome_occurrence'] = protein_df.apply(lambda row: [row['accession']+':'+str(row['consensus_epitopes_all'][i])+':'+str(row['core_epitopes_start_all'][i])+'-'+str(row['core_epitopes_end_all'][i]) for i in range(len(row['core_epitopes_start_all']))],axis=1)
    return protein_df


def reorder_peptides(row: pd.Series, intensity_column: str) -> pd.Series:
    """Reorder the peptides mapped to a protein by their position.
    
    Args:
        row: A row of a pandas dataframe containing per row one protein, 
            peptides mapped to the protein, start and end position of the peptide in the protein and the intensity of the peptide.
        intensity_column: The header of the column containing the intensities
            of the peptides.
    Returns:
        A reordered version of the input row, where the start positions are sorted in ascending order, the indices of the other columns are reordered in the same pattern. 
    """
    if intensity_column:
        lists = list(zip(row['start'], row['end'], row['sequence'], row['intensity']))
        sorted_lists = sorted(lists, key=lambda x: int(x[0]))
        starts, ends, sequences, intensities = zip(*sorted_lists)
        return list(starts), list(ends), list(sequences), list(intensities)
    else:
        lists = list(zip(row['start'], row['end'], row['sequence']))
        sorted_lists = sorted(lists, key=lambda x: int(x[0]))
        starts, ends, sequences = zip(*sorted_lists)
        return list(starts), list(ends), list(sequences)


def compute_consensus_epitopes(protein_df: pd.DataFrame, min_overlap: int, max_step_size: int, min_epi_len: int, intensity_column: float, mod_pattern: str) -> pd.DataFrame:
    """ Compute the core and whole sequence of all consensus epitope groups. 
    
    Args:
        protein_df: A pandas dataframe containing one protein per row.
        min_overlap: An integer of the minimal overlap between two epitopes     
            to be grouped to the same consensus epitope.
        max_step_size: An integer of the maximal distance between the start 
            position of two epitopes to be grouped to the same consensus 
            epitope.
        min_epi_len: An integer of the minimal length of a consensus epitope.
        intensity_column: The header of the column containing the intensities
            of the peptides.
        mod_pattern: A comma separated string with delimiters for peptide
            modifications

    Returns:
        The protein_df containing for each protein the core and whole sequence of each of its consensus epitope groups.
    """
    if intensity_column:
        protein_df[['start', 'end', 'sequence', 'intensity']] = protein_df.apply(lambda row: pd.Series(reorder_peptides(row, intensity_column)), axis=1)
    else:
        protein_df[['start', 'end', 'sequence']] = protein_df.apply(lambda row: pd.Series(reorder_peptides(row, intensity_column)), axis=1)
    # group peptides 
    protein_df = group_peptides(protein_df, min_overlap, max_step_size, intensity_column)
    protein_df = gen_landscape(protein_df,mod_pattern)
    protein_df = get_consensus_epitopes(protein_df, min_epi_len)
    return protein_df