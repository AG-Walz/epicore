"""
Reads in the evidence file and reference profile, computes the peptides positions in the proteome and links a protein accession with all peptides associated with the protein. 
"""

import pandas as pd 
import re
import numpy as np 
from Bio import SeqIO
import os 
import itertools

import logging
logger = logging.getLogger(__name__)


def read_id_output(id_output: str, seq_column: str, protacc_column: str, intensity_column: str, 
                   start_column: str, end_column: str, delimiter: str) -> pd.DataFrame:
    """Read in the evidence file.

    Args:
        id_output: The string of the path to the evidence file.
        seq_column: The string of the header of the column containing 
            peptide sequence information in the evidence file.
        protacc_column: The string of the header of the column containing 
            protein accession information in the evidence file.
        intensity_column: The string of the header of the column containing 
            intensity information in the evidence file.
        start_column: The string of the header of the column containing the 
            start positions of peptides in proteins.
        end_column: The string of the header of the column containing the end 
            position of peptides in proteins.
        delimiter: The delimiter that separates multiple entries in one column 
            in the evidence file.

    Returns:
        A pandas dataframe containing the columns sequence, protein accession, 
        peptides intensity and optional the columns start and end of the input 
        evidence file.

    Raises:
        Exception: If the file type of the provided evidence file is not 
            supported.
        Exception: If the provided column does not exist in the provided 
            evidence file.
        Exception: If a mandatory column header is not provided.
    """

    # determine the file type
    ext = os.path.splitext(id_output)[1]
    if ext == '.csv':
        peptides_df = pd.read_csv(id_output, delimiter=',')
    elif ext == '.tsv':
        peptides_df = pd.read_csv(id_output, delimiter='\t')
    elif ext == '.xlsx':
        peptides_df = pd.read_excel(id_output)
    else:
        raise Exception('The file type of your evidence file is not supported. \
                        Please use an evidence file that has one of the \
                        following file types: csv, tsv, xlsx')

    # check that the mandatory headers are provided and all provided column 
    # headers are part of the evidence file
    if seq_column not in peptides_df.columns:
        if seq_column:
            raise Exception('The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.'.
                            format(seq_column))
        else:
            raise Exception('The header for the column containing the peptide \
                            sequences in the evidence file is mandatory. \
                            Please provide the correct header name.')
    if protacc_column not in peptides_df.columns:
        if protacc_column:
            raise Exception('The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.'.
                            format(protacc_column))
        else:
            raise Exception('The header for the column containing the protein \
                            accessions in the evidence file is mandatory. \
                            Please provide the correct header name.')
    if (intensity_column not in peptides_df.columns) and intensity_column:
        raise Exception('The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.'.
                            format(intensity_column))
    if (start_column not in peptides_df.columns) and start_column:
        raise Exception('The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.'.
                            format(start_column))
    if (end_column not in peptides_df.columns) and end_column:
        raise Exception('The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.'.
                            format(end_column))

    # read in evidence file (start, stop position and intensity only when 
    # provided)
    if start_column and end_column:
        if intensity_column:
            peptides_df = peptides_df[[protacc_column,seq_column,       
                                       intensity_column, start_column, 
                                       end_column]]  
        else:
            peptides_df = peptides_df[[protacc_column,seq_column, start_column, 
                                       end_column]]
        peptides_df = peptides_df.astype(str)
        # split if peptide occurs multiple times in proteome 
        peptides_df[start_column] = peptides_df[start_column].apply(
            lambda x: list(x.split(delimiter)))
        peptides_df[end_column] = peptides_df[end_column].apply(
            lambda x: list(x.split(delimiter)))
        peptides_df[protacc_column] = peptides_df[protacc_column].apply(
            lambda x: list(x.split(delimiter)))
    else:
        if intensity_column:
            peptides_df = peptides_df[[protacc_column,seq_column, intensity_column]]
        else:
            peptides_df = peptides_df[[protacc_column,seq_column]]
        # split accessions if peptide occurs multiple times in proteome
        peptides_df[protacc_column] = peptides_df[protacc_column].apply(lambda x: list(set(x.split(delimiter))))

    return peptides_df


def proteome_to_dict(proteome: str) -> dict[str,str]:
    """Read reference proteome into dictionary.

    Args:
        proteome: The string of the path to the reference proteome.

    Returns: 
        The reference proteome as a dictionary.
    """
    proteome_dict = {}
    proteome = SeqIO.parse(open(proteome),'fasta')
    for protein in proteome:
        proteome_dict[protein.id] = str(protein.seq)
    return proteome_dict



def compute_pep_pos(peptide: str, accession: str, proteome_dict: dict[str,str]) -> tuple[list[int],list[int]]:
    """Compute the position of a peptide in a proteome.

    Args:
        peptide: The string of a peptide sequence.
        accession: The string of a protein accession.
        proteome_dict: A dictionary containing the reference proteome.

    Returns:
        Returns a tuple (start,end), where start is a list of all start positions of the peptide in the protein and end is a list of all end positions of the peptide in the protein.
    """   
    
    prot_seq = proteome_dict[accession]

    # regular expression for all occurrences of peptide in protein (also overlapping ones) 
    peptide_search = '(?=' + peptide + ')' 
    
    # get all occurrences of the peptide in the protein 
    pos_list = [[match.start(),match.start()+len(peptide)-1] for match in re.finditer(peptide_search,prot_seq)]
    if len(pos_list) == 0:
        raise Exception('The peptide {} does not occur in the protein with accession "{}"! Please use the proteome that was used for the identification of the peptides.'.format(peptide, accession))
    
    pos_list = np.array(pos_list)
    start = pos_list[:,0].astype(str).tolist()
    end = pos_list[:,1].astype(str).tolist()

    return start, end


def group_repetitive(starts: list[int], ends: list[int], peptide: str, accession:str)->tuple[list[int],list[int]]:
    """Group peptide occurrences that belong to the same repetitive region.

    Args: 
        peptide: The string of a peptide sequence.
        accession: The string of a protein accession.
        starts: A list of all start positions of the peptide in the protein.
        ends: A list of all end positions of the peptide in the protein. 

    Returns:
        Returns a tuple (starts,ends), where starts is a list containing the
        updated start positions and ends is a list of updated end positions. 
        These contain from the input start and end positions all start and end 
        positions that are not part of a repetitive region. For the start and 
        end positions that are part of repetitive regions the lowest start 
        position and highest end position is kept for each repetitive region. 
    """
    updated_starts = []
    updated_ends = []
    # add the first occurrences start positions to the start positions
    updated_starts.append(starts[0])
    for pep_pos in range(len(starts)-1):
        # two start positions are not part of one repetitive region if the next start position is higher than the current end position 
        if starts[pep_pos + 1] > ends[pep_pos]:
            updated_starts.append(starts[pep_pos + 1])
            updated_ends.append(ends[pep_pos])
    # add the last occurrences end position to the end positions
    updated_ends.append(ends[-1])
    if len(updated_starts) < len(starts): 
        return updated_starts, updated_ends
    else:
        return starts, ends


def get_start_end(row: pd.Series, peptide: str) -> tuple[list[int],list[int],list[float]]:
    """Gets the starts, ends and intensity of a peptide.
    
    Args:
        row: A row of a pandas dataframe containing per row one protein, 
            peptides mapped to the protein with their start and end position 
            and their intensity. 
        peptide: The string of a peptide sequence.

    Returns: 
        Returns a triple (starts, ends, intensity), where starts are the start 
        positions of the peptide in the protein, ends are the end positions of 
        the peptide in the protein and intensity is the intensity of the 
        peptide in the protein.
    """
    #TODO add correct headers and add intensity
    all_starts = row['start']
    all_ends = row['end']
    #intensities = row['intensities']
    peptides = row['sequence']
    starts = []
    ends = []
    for n_pep in range(len(peptides)):
        if peptide == peptides[n_pep]:
            starts.append(all_starts[n_pep])
            ends.append(all_ends[n_pep])
    return starts, ends
            


def prot_pep_link(peptides_df: pd.DataFrame, seq_column: str, protacc_column: str, intensity_column: str, start_column: str, end_column: str, proteome_dict: dict[str,str], mod_pattern:str) -> pd.DataFrame:
    """Converts a dataframe from one peptide per row to one protein per row.
    
    Args:
        peptides_df: A pandas dataframe containing the columns peptide sequence,
            protein accession, intensity, start position and end position.
        seq_column: The string of the header of the column containing 
            peptide sequence information in the evidence file.
        protacc_column: The string of the header of the column containing 
            protein accession information in the evidence file.
        intensity_column: The string of the header of the column containing 
            intensity information in the evidence file.
        start_column: The string of the header of the column containing the 
            start positions of peptides in proteins.
        end_column: The string of the header of the column containing the end 
            position of peptides in proteins.
        proteome_dict: A dictionary containing the reference proteome.
        mod_pattern: A comma separated string with delimiters for peptide
            modifications
    
    Returns:
        A pandas dataframe containing one protein per row and all peptides 
        mapped to that protein in the peptides_df, with their start position, 
        end position and intensity.

    Raises:
        Exception: If a peptide does not occur in the protein to which it is 
            mapped.
    """
    if not start_column and not end_column:
        # if the start and end positions of the peptides is not defined in the input evidence file 
        

        if intensity_column:
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'intensity'])
        else:
            proteins = pd.DataFrame(columns=['accession', 'sequence'])

        # create a row for each protein, that contains all peptides matched to the protein in the input evidence file
        for _, peptide in peptides_df.iterrows():
            # get all proteins matched to the peptide
            prot_accessions = peptide[protacc_column]
            for i, prot_accession in enumerate(prot_accessions):
                if prot_accession in proteins['accession'].values:
                    # update row for accessions seen before
                    proteins.loc[proteins['accession'] == prot_accession, 'sequence'] = proteins.loc[proteins['accession'] == prot_accession, 'sequence'].apply(lambda x: x + [peptide[seq_column]])
                    if intensity_column:
                        proteins.loc[proteins['accession'] == prot_accession, 'intensity'] = proteins.loc[proteins['accession'] == prot_accession, 'intensity'].apply(lambda x: x + [peptide[intensity_column]])

                else:
                    # create new row for accessions not seen before 
                    if intensity_column:
                        protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]], 'intensity':[peptide[intensity_column]]}
                    else:
                        protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]]}
                    proteins.loc[len(proteins)] = protein_entry

        # add the columns start and end
        proteins['start'] = [[] for _ in range(len(proteins))]
        proteins['end'] = [[] for _ in range(len(proteins))]

        # compute the peptide positions for each peptide of each protein
        for p, protein in proteins.iterrows():

            # get peptide sequence, protein accession and the measured intensity of the current row
            peptides = protein['sequence']
            accession = protein['accession']
            if intensity_column:
                intensities = protein['intensity']

            # starts and ends keep track of all occurrences of all peptides in peptides associated with the protein accession accession
            starts = []
            ends = []

            # updated_peps contains at each index the peptide associated with the start and end position in starts and ends at that index 
            # updated_intens contains at each index the intensity associated with the peptide at that index in updated_intens
            updated_peps = []
            if intensity_column:
                updated_intens = []

            for n_p, peptide in enumerate(peptides):
                # add each peptide and intensity associated with the accession to updated_peps and updated_intens
                updated_peps.append(peptide)
                if intensity_column:
                    updated_intens.append(intensities[n_p])
                
                # remove modifications in the sequence for the matching step (modifications are separated from the peptide by (), [] or by the user defined mod_pattern)
                pattern = r'\(.*?\)'
                peptide = re.sub(pattern,"",peptide)
                pattern = r'\[.*?\]'
                peptide = re.sub(pattern,"",peptide)
                if mod_pattern:
                    pattern = re.escape(mod_pattern.split(',')[0]) + r'.*?' + re.escape(mod_pattern.split(',')[1])
                    peptide = re.sub(pattern,"",peptide)
                pep_start, pep_end = compute_pep_pos(peptide, accession, proteome_dict)

                if len(pep_start) == 0:
                    raise Exception('The peptide {} does not occur in the protein with accession {} in the proteome you specified, but your input file provides evidence for that! Please use the proteome that was used for the identification of the peptides.'.format(peptide, prot_accession))
                if len(pep_start) > 1:
                    pep_start, pep_end = group_repetitive(pep_start,pep_end, peptide, accession)
                   
                if len(pep_start) > 1:
                    # if a peptide occurs multiple times in a protein
                    for _ in pep_start[:-1]:
                        updated_peps.append(peptide)
                        if intensity_column:
                            updated_intens.append(intensities[n_p])
                
                # collect all start and end positions of the peptide in the protein
                starts.extend(pep_start)
                ends.extend(pep_end)
            
            proteins.at[p, 'start'] = starts
            proteins.at[p, 'end'] = ends
            proteins.at[p, 'sequence'] = updated_peps
            if intensity_column:
                proteins.at[p,'intensity'] = updated_intens
    else:
        if intensity_column:
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'intensity', 'start','end'])
        else:
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'start','end'])

        for _, peptide in peptides_df.iterrows():
            # get all proteins associated with the peptide
            prot_accessions = peptide[protacc_column]
            for i, prot_accession in enumerate(prot_accessions):

                if prot_accession in proteins['accession'].values:
                    # append sequence, intensity, start and end information for accessions seen before
                    proteins.loc[proteins['accession'] == prot_accession, 'sequence'] = proteins.loc[proteins['accession'] == prot_accession, 'sequence'].apply(lambda x: x + [peptide[seq_column]])
                    if intensity_column:
                        proteins.loc[proteins['accession'] == prot_accession, 'intensity'] = proteins.loc[proteins['accession'] == prot_accession, 'intensity'].apply(lambda x: x + [peptide[intensity_column]])
                    proteins.loc[proteins['accession'] == prot_accession, 'start'] = proteins.loc[proteins['accession'] == prot_accession, 'start'].apply(lambda x: x + [peptide[start_column][i]])
                    proteins.loc[proteins['accession'] == prot_accession, 'end'] = proteins.loc[proteins['accession'] == prot_accession, 'end'].apply(lambda x: x + [peptide[end_column][i]])
                    
                else:
                    # create new row for accessions not seen before
                    if intensity_column:
                        protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]], 'intensity':[peptide[intensity_column]], 'start':[peptide[start_column][i]], 'end':[peptide[end_column][i]]}
                    else:
                        protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]], 'start':[peptide[start_column][i]], 'end':[peptide[end_column][i]]}
                    proteins.loc[len(proteins)] = protein_entry

        
        for p, protein in proteins.iterrows():

            # get peptide sequence, protein accession and the measured intensity of the current row
            peptides = protein['sequence']
            accession = protein['accession']
            if intensity_column:
                intensities = protein['intensity']
            starts = []
            ends = []
            
            # updated_peps contains at each index the peptide associated with the start and end position in starts and ends at that index 
            # updated_intens contains at each index the intensity associated with the peptide at that index in updated_intens
            updated_peps = []
            if intensity_column:
                updated_intens = []

            pep_pos_seen = []
            for n_p, peptide in enumerate(peptides):
                if peptide in pep_pos_seen:
                    continue
                else: 
                    pep_pos_seen.append(peptide)

                # add each peptide and intensity associated with the accession to updated_peps and updated_intens
                updated_peps.append(peptide)
                if intensity_column:
                    updated_intens.append(intensities[n_p])
                
                # get all starts, ends and the intensity of that peptide in that row
                pep_start, pep_end = get_start_end(protein, peptide)
            
                if len(pep_start) == 0:
                    raise Exception('The peptide {} does not occur in the protein with accession {} in the proteome you specified, but your input file provides evidence for that! Please use the proteome that was used for the identification of the peptides.'.format(peptide, prot_accession))
                if len(pep_start) > 1:
                    pep_start, pep_end = group_repetitive(pep_start,pep_end, peptide, accession)
                    
                   
                if len(pep_start) > 1:
                    # if a peptide occurs multiple times in a protein
                    for _ in pep_start[:-1]:
                        updated_peps.append(peptide)
                        if intensity_column:
                            updated_intens.append(intensities[n_p])
                
                # collect all start and end positions of the peptide in the protein
                starts.extend(pep_start)
                ends.extend(pep_end)
            

            proteins.at[p, 'start'] = starts
            proteins.at[p, 'end'] = ends
            proteins.at[p, 'sequence'] = updated_peps
            if intensity_column:
                proteins.at[p,'intensity'] = updated_intens
    return proteins


def parse_input(evidence_file: str, seq_column: str, protacc_column: str, intensity_column: str, start_column: str, end_column: str, delimiter: str, proteome_dict: dict[str,str], mod_pattern: str) -> pd.DataFrame:
    """Parse the evidence file.
    
    Args:
        evidence_file: The string of the path to the evidence file.
        seq_column: The string of the header of the column containing 
            peptide sequence information in the evidence file.
        protacc_column: The string of the header of the column containing 
            protein accession information in the evidence file.
        intensity_column: The string of the header of the column containing 
            intensity information in the evidence file.
        start_column: The string of the header of the column containing the 
            start positions of peptides in proteins.
        end_column: The string of the header of the column containing the end 
            position of peptides in proteins.
        delimiter: The delimiter that separates multiple entries in one column 
            in the evidence file.
        proteome_dict: A dictionary containing the reference proteome.
        mod_pattern: A comma separated string with delimiters for peptide
            modifications
    
    Returns:
        A tuple of a pandas dataframe and an integer. The dataframe contains 
        one protein per row and all peptides mapped to that protein in the 
        peptides_df, with their start position, end position and intensity.
        The integer is the number of peptides that are removed due to their 
        accession not appearing in the proteome file. 
    """
    peptides_df = read_id_output(evidence_file, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter)

    # get peptides/proteins with protein accessions that do not appear in the proteome
    peptides = peptides_df.apply(lambda row: [prot for prot in row[protacc_column] if prot not in proteome_dict.keys()], axis=1).values
    n_removed_proteins = set(itertools.chain.from_iterable(peptides))

    # remove peptides with protein accessions that do not appear in the proteome 
    if start_column and end_column:
        peptides_df[start_column] = peptides_df.apply(lambda row: [start for start, prot in zip(row[start_column],row[protacc_column]) if prot in proteome_dict.keys()], axis=1)
        peptides_df[end_column] = peptides_df.apply(lambda row: [end for end, prot in zip(row[end_column],row[protacc_column]) if prot in proteome_dict.keys()], axis=1)
    peptides_df[protacc_column] = peptides_df.apply(lambda row: [prot for prot in row[protacc_column] if prot in proteome_dict.keys()], axis=1)

    logger.info(f'Peptides mapped to the following {len(n_removed_proteins)} proteins were removed since the proteins do not appear in the proteome fasta file: {n_removed_proteins}.')
    protein_df = prot_pep_link(peptides_df, seq_column, protacc_column, intensity_column, start_column, end_column, proteome_dict, mod_pattern)
    return protein_df, len(n_removed_proteins)
