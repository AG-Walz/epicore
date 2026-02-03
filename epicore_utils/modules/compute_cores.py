"""
Computes core epitopes, by grouping overlapping peptides, together, building a landscape for each group and identifying plateaus with a defined minimal length in each landscape. 
"""

import numpy as np
import re
import pandas as pd
from multiprocessing import get_context, cpu_count


def group_peptides_protein(row: list, min_overlap: int, max_step_size: int, intensity_column: str, total_intens: float, strict: bool) -> pd.DataFrame:
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
        total_intens: The total intensity of the evidence file.
        strict: Indicates if strict version should be run.

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

    start_pos =  row['start']
    end_pos = row['end']
    sequences = row['sequence']
    samples = row['sample']
    conditions = row['condition']
    if intensity_column:
        intensity = row['intensity']
    peptide_indices = row['peptide_index']

    grouped_peptides_start = []
    grouped_peptides_end = []
    grouped_peptides_sequence = []
    grouped_peptides_sample = []
    grouped_peptides_condition = []
    grouped_peptides_max_end = 0
    grouped_peptides_index = []
    if intensity_column:
        grouped_peptides_intensity = []
        core_intensity = 0
    n_jumps = 0
    mapping = []

    for i in range(len(start_pos)-1):

        grouped_peptides_start.append(int(start_pos[i]))
        grouped_peptides_end.append(int(end_pos[i]))
        grouped_peptides_sequence.append(sequences[i])
        grouped_peptides_sample.append(samples[i])
        grouped_peptides_condition.append(conditions[i])
        grouped_peptides_max_end = max(grouped_peptides_max_end, int(end_pos[i]))
        grouped_peptides_index.append(peptide_indices[i])
        if intensity_column:
            grouped_peptides_intensity.append(intensity[i])

        step_size = int(start_pos[i+1]) - int(start_pos[i])
        mapping.append(n_jumps)
        if intensity_column:
            core_intensity += float(intensity[i])

        #first_start = grouped_peptides_start[0]
        first_end = grouped_peptides_end[0]
        group_overlap = max(min(int(end_pos[i+1])-int(start_pos[i+1])+1, first_end-int(start_pos[i+1])+1),0) 

        overlap = int(end_pos[i]) - int(start_pos[i+1]) +1
        group_overlap = min(grouped_peptides_end) - int(start_pos[i+1]) +1

        if strict:
            condition_group = (((step_size >= max_step_size) and (overlap < min_overlap)) or (overlap < min_overlap) or (group_overlap < min_overlap)) and (grouped_peptides_max_end<int(end_pos[i+1]))
        else:
            condition_group = ((step_size >= max_step_size) and (overlap < min_overlap))

        # create new peptide group after each jump
        if step_size != 0:

            if condition_group:
                row['grouped_peptides_start'].append(grouped_peptides_start)
                row['grouped_peptides_end'].append(grouped_peptides_end)
                row['grouped_peptides_sequence'].append(grouped_peptides_sequence)
                row['grouped_peptides_sample'].append(grouped_peptides_sample)
                row['grouped_peptides_condition'].append(grouped_peptides_condition)
                row['peptide_indices'].append(grouped_peptides_index)
                if intensity_column:
                    row['grouped_peptides_intensity'].append(grouped_peptides_intensity)
                    row['core_epitopes_intensity'].append(core_intensity)                   
                    
                n_jumps += 1
                grouped_peptides_end = []
                grouped_peptides_start = []
                grouped_peptides_sequence = []
                grouped_peptides_sample = []
                grouped_peptides_condition = []
                grouped_peptides_intensity = []
                grouped_peptides_index = []
                grouped_peptides_max_end = 0
                if intensity_column:
                    core_intensity = 0


    # special case for last peptide match of protein
    if len(grouped_peptides_end) == 0:

        row['grouped_peptides_start'].append([int(start_pos[-1])])
        row['grouped_peptides_end'].append([int(end_pos[-1])])
        row['grouped_peptides_sequence'].append([sequences[-1]])
        row['grouped_peptides_sample'].append([samples[-1]])
        row['grouped_peptides_condition'].append([conditions[-1]])
        row['peptide_indices'].append([peptide_indices[-1]])
        if intensity_column:
            row['grouped_peptides_intensity'].append([intensity[-1]])
            row['core_epitopes_intensity'].append(peptide_indices[-1])
        mapping.append(n_jumps)
    
    else:

        grouped_peptides_start.append(int(start_pos[-1]))
        grouped_peptides_end.append(int(end_pos[-1]))
        grouped_peptides_sequence.append(sequences[-1])
        grouped_peptides_sample.append(samples[-1])
        grouped_peptides_condition.append(conditions[-1])
        grouped_peptides_index.append(peptide_indices[-1])
        if intensity_column:
            grouped_peptides_intensity.append(intensity[-1])
            core_intensity += float(intensity[-1])
            row['grouped_peptides_intensity'].append(grouped_peptides_intensity)
            row['core_epitopes_intensity'].append(core_intensity)
            

        row['grouped_peptides_start'].append(grouped_peptides_start)
        row['grouped_peptides_end'].append(grouped_peptides_end)
        row['grouped_peptides_sequence'].append(grouped_peptides_sequence)
        row['grouped_peptides_sample'].append(grouped_peptides_sample)
        row['grouped_peptides_condition'].append(grouped_peptides_condition)
        row['peptide_indices'].append(grouped_peptides_index)
        mapping.append(n_jumps)

    row['sequence_group_mapping'] = mapping

    # update start and end position (add start and end position for multiple occurrences)
    row['start'] = [f'{start}' for group in row['grouped_peptides_start'] for start in group]
    row['end'] = [f'{end}' for group in row['grouped_peptides_end'] for end in group]
    row['peptide_index'] = [f'{index}' for group in row['peptide_indices'] for index in group]

    return row

def group_peptides(protein_df: pd.DataFrame, min_overlap: int, max_step_size: int, intensity_column: str, total_intens: float, strict: bool) -> pd.DataFrame:
    
    # start, end, sequence and intensity of peptides of one group grouped together
    protein_df['grouped_peptides_start'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_end'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_sequence'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_sample'] = [[] for _ in range(len(protein_df))]
    protein_df['grouped_peptides_condition'] = [[] for _ in range(len(protein_df))]
    if intensity_column:
        protein_df['grouped_peptides_intensity'] = [[] for _ in range(len(protein_df))]
        # for each peptide group the total and relative intensity of the entire group
        protein_df['core_epitopes_intensity'] = [[] for _ in range(len(protein_df))]
        protein_df['relative_core_intensity'] = [[] for _ in range(len(protein_df))]
    # for each peptide the index of its group
    protein_df['sequence_group_mapping'] = [[] for _ in range(len(protein_df))]
    protein_df['peptide_indices'] = [[] for _ in range(len(protein_df))]
    
    protein_df = protein_df.apply(lambda row: group_peptides_protein(row, min_overlap, max_step_size, intensity_column, total_intens, strict), axis=1)
    return protein_df


def pep_landscape(row: pd.Series) -> np.array:
    """Generate landscape of a peptide in a group.

    Args:
        row: A row of a pandas dataframe containing per row one protein, 
             peptides mapped to the protein with their start and end position 
             and their intensity. 

    Returns:
        A list containing the landscape of a protein in a given peptide group.
    """
    zeros = np.zeros([row['end_max']-row['start_min']+1])
    zeros[row['grouped_peptides_start']-row['start_min']:(row['grouped_peptides_end']-row['start_min']+1)] = 1
    return zeros

def landscape_chunk(chunk_df):
    """Generate landscapes of the peptides in chunk_df.

    Args:
        chunk_df: A subset of the protein DataFrame.

    Returns:
        The input DataFrame with the additional column pep_landscape.
    """
    chunk_df['pep_landscape'] = chunk_df.apply(lambda row: pep_landscape(row), axis=1)
    return chunk_df

def comp_landscape(protein_df: pd.DataFrame, proteome_dict: dict[str,str]) -> pd.DataFrame:
    """Compute the landscape of all consensus epitope groups.
    
    Args:
        protein_df: A pandas dataframe containing one protein per row.
        proteome_dict: A dictionary containing the reference proteome.

    Returns:
        The protein_df with two additional columns (landscape, whole_epitopes).
        The column landscape contains the landscapes of the consensus epitopes 
        group. The height of the landscape of a consensus epitope group at a 
        position is the number of peptides that contain that position. The 
        column whole_epitopes contains the entire sequence of the consensus 
        epitope group.
    """
    # get entire sequence of each group
    protein_df['start_min'] = [min(x) for x in protein_df['grouped_peptides_start']]
    protein_df['end_max'] = [max(x) for x in protein_df['grouped_peptides_end']]
    protein_df['whole_epitopes'] = protein_df['accession'].map(proteome_dict)    
    protein_df['whole_epitopes'] = [whole_epitopes[start:end+1] for start, end, whole_epitopes in zip(protein_df['start_min'],protein_df['end_max'],protein_df['whole_epitopes'])]#protein_df['whole_epitopes'].str.slice(start=protein_df['start_min'],stop=protein_df['end_max']+1)
    protein_df = protein_df.assign(whole_epitopes_all=protein_df.whole_epitopes)

    protein_df['group'] = range(len(protein_df))
    protein_df = protein_df.explode(['grouped_peptides_start', 'grouped_peptides_end', 'grouped_peptides_sample', 'grouped_peptides_condition', 'grouped_peptides_sequence'])

    n_parallel = max(1, cpu_count()-5)
    # compute the landscape for each peptide in the group
    blocksize = max(len(protein_df) // n_parallel,1)
    with get_context('spawn').Pool(min(n_parallel,blocksize)) as pool:
        chunk_dfs = pool.starmap(landscape_chunk,[(protein_df.iloc[chunk*blocksize:(chunk+1)*blocksize if chunk < min(n_parallel,blocksize) -1 else len(protein_df)],) for chunk in range(min(n_parallel,blocksize))])
    protein_df = pd.concat(chunk_dfs)

    cols = protein_df.columns.values
    aggr = {'whole_epitopes_all':lambda x: list(x), 'landscape': 'first', 'grouped_peptides_start': lambda x: list(x), 'grouped_peptides_end': lambda x: list(x), 'grouped_peptides_sample': lambda x: list(x), 'grouped_peptides_condition': lambda x: list(x), 'grouped_peptides_sequence': lambda x: list(x)}
    for col in cols:
        if col not in aggr.keys():
            aggr[col] = 'first'

    comb_landscapes = protein_df.groupby('group').agg({'pep_landscape': lambda column: np.sum(tuple(column), axis=0)}).reset_index()
    comb_landscapes = comb_landscapes.rename(columns={'pep_landscape':'landscape'})

    # compute the landscape of the entire group
    protein_df = pd.merge(protein_df, comb_landscapes, on='group')
    protein_df['landscape'] = protein_df['landscape'].apply(lambda landscape: landscape.tolist())
    protein_df = protein_df.groupby('group').agg(aggr)
    aggr['landscape'] = lambda x: list(x)
    aggr['whole_epitopes'] = lambda x: list(x)
    aggr['whole_epitopes'] = lambda x: list(x)
    del aggr['accession']
    aggr['start_min'] = lambda x: list(x)

    # combine all groups of one protein
    protein_df = protein_df.groupby('accession').agg(aggr)
    protein_df['whole_epitopes_all'] = protein_df['whole_epitopes_all'].apply(lambda x: np.hstack(x).tolist())
    protein_df = protein_df.drop(['group','pep_landscape'], axis=1)
    protein_df = protein_df.reset_index()

    return protein_df


def find_minima(landscape: np.array) -> np.array:
    """Compute local minima of a landscape.
    
    Args:
        landscape: A list containing the landscape of a peptide group.

    Returns:
        A list of the positions where the landscape has a minimum
    """
    # find minima
    landscape_f = np.roll(landscape, -1) > landscape
    landscape_b = np.roll(landscape, 1) > landscape

    # ensure minima is not at the edge
    edge = np.full(len(landscape), True)
    edge[0] = False
    edge[-1] = False
    landscape = landscape_f & landscape_b & edge
    return landscape


def new_groups(starts: list[int], ends: list[int], samples: list[str], conditions: list[str], splits: list[int], sequences: list[str]) -> list:
    """Generate new peptide groups based on the landscape minima. 

    Args:
        starts: A list containing all start positions of the current peptide group. 
        ends: A list containing all end positions of the current peptide group. 
        samples: A list containing all samples of the current peptide group. 
        conditions: A list containing all conditions of the current peptide group. 
        splits: A list containing all minima of the landscape of the current peptide group. 
        sequences: A list containing all sequences of the current peptide group.
    
    Returns:
        A list containing the new peptide groups.
    """
    start_groups = []
    end_groups = []
    sample_groups = []
    condition_groups = []
    sequence_groups = []
    start_list = []
    end_list = []
    sample_list = []
    condition_list = []
    sequence_list = []

    # if the landscape has no minima
    if len(splits) == 0:
        return [[starts], [ends], [sequences], [samples], [conditions]]

    
    for start, end, sample, condition, sequence in zip(starts, ends, samples, conditions, sequences):
        
        if len(splits) == 0 or splits[0] > start:
            start_list.append(start)
            end_list.append(end)
            sample_list.append(sample)
            condition_list.append(condition)
            sequence_list.append(sequence)

        elif splits[0] <= start: # TODO check if split position is correct here (right index and < or <=)
            # start new peptide group after minimum
            start_groups.append(start_list)
            end_groups.append(end_list)
            sample_groups.append(sample_list)
            condition_groups.append(condition_list)
            sequence_groups.append(sequence_list)
            sequence_list = [sequence]
            start_list = [start]
            end_list = [end]
            sample_list = [sample]
            condition_list = [condition]
            splits = splits[1:]

    start_groups.append(start_list)
    end_groups.append(end_list)
    sample_groups.append(sample_list)
    condition_groups.append(condition_list)
    sequence_groups.append(sequence_list)

    return [start_groups, end_groups, sequence_groups, sample_groups, condition_groups]


def group_refinement(protein_df: pd.DataFrame, proteome_dict: dict[str,str]):
    """Split peptide groups at landscape minima.

    Args: 
        protein_df: A pandas dataframe containing one protein per row.
        proteome_dict: A dictionary containing the reference proteome.

    Returns:
        The pandas dataframe with refined peptide groups.
    """

    # convert protein_df so each peptide group is in one row
    protein_df = protein_df.explode(['grouped_peptides_start', 'grouped_peptides_end', 'grouped_peptides_sample', 'grouped_peptides_condition', 'grouped_peptides_sequence', 'landscape', 'start_min'])
    
    # find landscape minima
    protein_df['split_rows'] = protein_df['landscape'].apply(lambda landscape: find_minima(landscape))
    protein_df['split_position'] = protein_df.apply(lambda row: np.where(row['split_rows'])[0]+row['start_min'], axis=1)

    # refine peptide groups
    protein_df[['new_groups_start', 'new_groups_end', 'new_group_sequences', 'new_groups_sample', 'new_groups_condition']] = protein_df.apply(lambda row: new_groups(row['grouped_peptides_start'], row['grouped_peptides_end'], row['grouped_peptides_sample'], row['grouped_peptides_condition'], row['split_position'], row['grouped_peptides_sequence']), axis=1, result_type='expand')
    protein_df = protein_df.explode(['new_groups_start', 'new_groups_end', 'new_groups_sample', 'new_groups_condition', 'new_group_sequences'])
    protein_df = protein_df.drop(['grouped_peptides_start', 'grouped_peptides_end', 'grouped_peptides_sample', 'grouped_peptides_condition', 'grouped_peptides_sequence'], axis=1)
    protein_df = protein_df.rename(columns={'new_groups_start':'grouped_peptides_start', 'new_groups_end':'grouped_peptides_end', 'new_groups_sample':'grouped_peptides_sample','new_groups_condition':'grouped_peptides_condition', 'new_group_sequences': 'grouped_peptides_sequence'})
    protein_df = protein_df.drop(['landscape'], axis=1)

    # recompute group landscapes
    #protein_df = protein_df.group_by(['accession'])
    protein_df = comp_landscape(protein_df,proteome_dict)

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
                    pep_in_prot_start = ce_start_pos.item()
                    pep_in_prot_end = pep_in_prot_start + current_pep_length.item()

                    # get consensus epitopes
                    whole_epitope_wo_mod = protein_df.at[r,'whole_epitopes'][group]
                    for _ in row['grouped_peptides_sequence'][group]:
                        protein_df.at[r,'consensus_epitopes_all'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                        protein_df.at[r,'core_epitopes_start_all'].append((pep_in_prot_start+min(row['grouped_peptides_start'][group])))
                        protein_df.at[r,'core_epitopes_end_all'].append(pep_in_prot_end+min(row['grouped_peptides_start'][group]) - 1) 
                    protein_df.at[r,'consensus_epitopes'].append(whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end])
                    protein_df.at[r,'core_epitopes_start'].append(pep_in_prot_start+min(row['grouped_peptides_start'][group]))
                    protein_df.at[r,'core_epitopes_end'].append(pep_in_prot_end+min(row['grouped_peptides_start'][group]) - 1)
                    break
                
                # if no core with length > min_epi_length
                if total_count == total_counts[-1]:
                    pep_in_prot_start = ce_start_pos.item()
                    pep_in_prot_end = pep_in_prot_start + current_pep_length.item()
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
        lists = list(zip(row['start'], row['end'], row['sequence'], row['intensity'], row['peptide_index'], row['sample'], row['condition']))
        lists = sorted(lists, key=lambda x: x[5], reverse=True)
        lists = sorted(lists, key=lambda x: int(x[1]), reverse=True)
        lists = sorted(lists, key=lambda x: x[2], reverse=True)
        sorted_lists = sorted(lists, key=lambda x: int(x[0]))
        starts, ends, sequences, intensities, indices, samples, conditions = zip(*sorted_lists)
        return list(starts), list(ends), list(sequences), list(intensities), list(indices), list(samples), list(conditions)
    else:
        lists = list(zip(row['start'], row['end'], row['sequence'], row['peptide_index'], row['sample'], row['condition']))
        lists = sorted(lists, key=lambda x: x[4], reverse=True)
        lists = sorted(lists, key=lambda x: x[2], reverse=True)
        lists = sorted(lists, key=lambda x: int(x[1]), reverse=True)
        sorted_lists = sorted(lists, key=lambda x: int(x[0]))
        starts, ends, sequences, indices, samples, conditions = zip(*sorted_lists)
        return list(starts), list(ends), list(sequences), list(indices), list(samples), list(conditions)


def compute_consensus_epitopes(protein_df: pd.DataFrame, min_overlap: int, max_step_size: int, min_epi_len: int, intensity_column: float, mod_pattern: str, proteome_dict: dict[str,str], total_intens: float, strict: bool) -> pd.DataFrame:
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
        proteome_dict: A dictionary containing the reference proteome.
        total_intens: The total intensity of the evidence file.
        strict: A boolean indicating if the strict version should be run.

    Returns:
        The protein_df containing for each protein the core and whole sequence of each of its consensus epitope groups.
    """
    if intensity_column:
        protein_df[['start', 'end', 'sequence', 'intensity', 'peptide_index', 'sample', 'condition']] = protein_df.apply(lambda row: pd.Series(reorder_peptides(row, intensity_column)), axis=1)
    else:
        protein_df[['start', 'end', 'sequence', 'peptide_index', 'sample', 'condition']] = protein_df.apply(lambda row: pd.Series(reorder_peptides(row, intensity_column)), axis=1)
    # group peptides
    protein_df = group_peptides(protein_df, min_overlap, max_step_size, intensity_column, total_intens, strict)
    protein_df = protein_df.explode(['grouped_peptides_start', 'grouped_peptides_end', 'grouped_peptides_sample', 'grouped_peptides_sequence', 'grouped_peptides_condition'])
    protein_df = comp_landscape(protein_df, proteome_dict)
    if not strict:
        protein_df = protein_df[['accession','sequence','start','end','peptide_index','sample','condition','grouped_peptides_start','grouped_peptides_end','grouped_peptides_sequence','grouped_peptides_sample', 'grouped_peptides_condition','sequence_group_mapping','landscape','whole_epitopes','whole_epitopes_all', 'start_min']]
        protein_df = group_refinement(protein_df, proteome_dict)
    protein_df = protein_df[['accession','sequence','start','end','peptide_index','sample','condition','grouped_peptides_start','grouped_peptides_end','grouped_peptides_sequence','grouped_peptides_sample', 'grouped_peptides_condition','sequence_group_mapping','landscape','whole_epitopes','whole_epitopes_all']]
    protein_df = get_consensus_epitopes(protein_df, min_epi_len)
    return protein_df
