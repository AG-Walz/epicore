"""
Reads in the evidence file and reference profile, computes the peptides positions in the proteome and links a protein accession with all peptides associated with the protein. 
"""

import pandas as pd 
import re
import numpy as np 
from Bio import SeqIO
import os 
import itertools
import polars as pl

import logging
logger = logging.getLogger(__name__)


def read_id_output(id_output: str, seq_column: str, protacc_column: str, intensity_column: str, 
                   start_column: str, end_column: str, delimiter: str, sample_column: str, condition_column: str) -> pd.DataFrame:
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
        peptides_df = pl.read_csv(id_output, separator=',').with_row_count('peptide_index')
    elif ext == '.tsv':
        peptides_df = pl.read_csv(id_output, separator='\t').with_row_count('peptide_index')
    elif ext == '.xlsx':
        peptides_df = pl.read_excel(id_output)
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
            peptides_df = peptides_df.select(protacc_column, seq_column, intensity_column, 'peptide_index', start_column, end_column, sample_column, condition_column)
        else:
            peptides_df = peptides_df.select(protacc_column, seq_column, 'peptide_index', start_column, end_column, sample_column, condition_column)

        # split if peptide occurs multiple times in proteome
        peptides_df = peptides_df.with_columns(pl.col(start_column).cast(pl.String))
        peptides_df = peptides_df.with_columns(pl.col(end_column).cast(pl.String))
        peptides_df = peptides_df.with_columns((pl.col(start_column).str.split(delimiter).alias(start_column)))
        peptides_df = peptides_df.with_columns((pl.col(end_column).str.split(delimiter).alias(end_column)))
        peptides_df = peptides_df.with_columns((pl.col(protacc_column).str.split(delimiter).alias(protacc_column)))
    
    else:
        if intensity_column:
            peptides_df = peptides_df.select(protacc_column, seq_column, intensity_column, 'peptide_index', sample_column, condition_column)
        else:
            peptides_df = peptides_df.select(protacc_column, seq_column, 'peptide_index', sample_column, condition_column)
        
        # split accessions if peptide occurs multiple times in proteome
        peptides_df = peptides_df.with_columns((pl.col(protacc_column).str.split(delimiter)).alias(protacc_column))

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


def group_repetitive(starts: list[int], ends: list[int], peps: list[str], accs, idex, samples, conditions)->tuple[list[int],list[int]]:
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
    updated_pos = [] # TODO have to be sorted

    for start, end, pep, acc, idx, sample, condition in zip(starts,ends, peps, accs, idex, samples, conditions):
        current = -1
        updated_start = []
        updated_end = []
        updated_idx = []
        updated_conditions = []
        # add the first occurrences start positions to the start positions
        updated_start.append(str(start[0]))
        updated_peps = []
        updated_samples = []
        lists = list(zip(start, end, idx))
        lists = sorted(lists, key=lambda x: int(x[0]))
        start, end, idx = zip(*lists)
        group_ends = []
        for pep_pos in range(len(start)-1):

            # two start positions are not part of one repetitive region if the next start position is higher than the current end position 
            if int(start[pep_pos + 1]) > int(end[pep_pos]): # new group
                group_ends.append(end[pep_pos])
                # add max end of repetitive group
                for group_end in group_ends:
                    updated_end.append(str(max(group_ends)))
                group_ends = []
                updated_start.append(str(start[pep_pos + 1]))
                updated_idx.append(str(idx[pep_pos].to_list()[0]))
                updated_peps.append(pep)
                updated_samples.append(sample[0])
                updated_conditions.append(condition[0])
                current = pep_pos

            else: #  add min_start for repetitive group
                group_ends.append(end[pep_pos])
                updated_start.append(str(start[current + 1]))
                updated_idx.append(str(idx[current].to_list()[0]))
                updated_peps.append(pep)
                updated_samples.append(sample[0])
                updated_conditions.append(condition[0])

        # add the last occurrences end position to the end positions
        if len(group_ends) == 0:
            updated_end.append(str(end[-1]))
        else:
            group_ends.append(end[-1])
            for group_end in group_ends:
                updated_end.append(str(max(group_ends)))
            #updated_end = updated_end[:-1]

        updated_idx.append(str(idx[-1].to_list()[0]))
        updated_peps.append(pep)
        updated_samples.append(sample[0])
        updated_conditions.append(condition[0])

        # reduce each occurrence to one
        updated_df = pd.DataFrame({'start':updated_start, 'end':updated_end, 'peps':updated_peps, 'idx':updated_idx, 'sample':updated_samples, 'cond':updated_conditions})
        updated_df = updated_df.drop_duplicates()
        updated_start = ';'.join(updated_df['start'])
        updated_end = ';'.join(updated_df['end'])
        updated_peps = ';'.join(updated_df['peps'])
        updated_idx = ';'.join(updated_df['idx'])
        updated_samples = ';'.join(updated_df['sample'])
        updated_conditions = ';'.join(updated_df['cond'])

        updated_pos.append(f'{updated_start}|{updated_end}|{updated_idx}|{updated_peps}|{updated_samples}|{updated_conditions}')

    return updated_pos


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
            


def prot_pep_link(peptides_df: pd.DataFrame, seq_column: str, protacc_column: str, intensity_column: str, start_column: str, end_column: str, proteome_dict: dict[str,str], mod_pattern:str, delimiter, sample_column, condition_column) -> pd.DataFrame:
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
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'intensity', 'peptide_index'])
        else:
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'peptide_index'])

        # create a row for each protein, that contains all peptides matched to the protein in the input evidence file
        for index_peptide, peptide in peptides_df.iterrows():
            
            # get all proteins matched to the peptide
            prot_accessions = peptide[protacc_column]
            for i, prot_accession in enumerate(prot_accessions):
                if prot_accession in proteins['accession'].values:
                    # update row for accessions seen before
                    proteins.loc[proteins['accession'] == prot_accession, 'sequence'] = proteins.loc[proteins['accession'] == prot_accession, 'sequence'].apply(lambda x: x + [peptide[seq_column]])
                    proteins.loc[proteins['accession'] == prot_accession, 'peptide_index'] = proteins.loc[proteins['accession'] == prot_accession, 'peptide_index'].apply(lambda x: x + [index_peptide])
                    if intensity_column:
                        proteins.loc[proteins['accession'] == prot_accession, 'intensity'] = proteins.loc[proteins['accession'] == prot_accession, 'intensity'].apply(lambda x: x + [peptide[intensity_column]])

                else:
                    # create new row for accessions not seen before 
                    if intensity_column:
                        protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]], 'intensity':[peptide[intensity_column]], 'peptide_index':[index_peptide]}
                    else:
                        protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]], 'peptide_index':[index_peptide]}
                    proteins.loc[len(proteins)] = protein_entry

        # add the columns start and end
        proteins['start'] = [[] for _ in range(len(proteins))]
        proteins['end'] = [[] for _ in range(len(proteins))]

        # compute the peptide positions for each peptide of each protein
        for p, protein in proteins.iterrows():

            # get peptide sequence, protein accession and the measured intensity of the current row
            peptides = protein['sequence']
            accession = protein['accession']
            indices = protein['peptide_index']
            if intensity_column:
                intensities = protein['intensity']

            # starts and ends keeps track of all occurrences of all peptides in peptides associated with the protein accession
            starts = []
            ends = []

            # updated_peps contains at each index the peptide associated with the start and end position in starts and ends at that index 
            # updated_intens contains at each index the intensity associated with the peptide at that index in updated_intens
            updated_peps = []
            updated_index = []
            if intensity_column:
                updated_intens = []

            for n_p, peptide in enumerate(peptides):
                # add each peptide and intensity associated with the accession to updated_peps and updated_intens
                updated_peps.append(peptide)
                updated_index.append(indices[n_p])
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
                        updated_index.append(indices[n_p])
                        if intensity_column:
                            updated_intens.append(intensities[n_p])
                
                # collect all start and end positions of the peptide in the protein
                starts.extend(pep_start)
                ends.extend(pep_end)
            
            proteins.at[p, 'start'] = starts
            proteins.at[p, 'end'] = ends
            proteins.at[p, 'sequence'] = updated_peps
            proteins.at[p, 'peptide_index'] = updated_index
            if intensity_column:
                proteins.at[p,'intensity'] = updated_intens
    else:
        if intensity_column:
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'intensity', 'start','end', 'peptide_index'])
        else:
            proteins_df = peptides_df.explode(protacc_column, start_column, end_column)
            proteins_df = proteins_df.with_columns(pl.col(start_column).cast(pl.Int64))
            proteins_df = proteins_df.sort(start_column).group_by(protacc_column).agg(pl.col(seq_column), pl.col(start_column), pl.col(end_column), pl.col('peptide_index'), pl.col(sample_column), pl.col(condition_column))
            proteins_df = proteins_df.explode(start_column, end_column, seq_column, 'peptide_index', sample_column, condition_column).group_by(seq_column, protacc_column, sample_column, condition_column).agg(pl.col(start_column), pl.col(end_column), pl.col('peptide_index'))
            proteins_df = proteins_df.with_columns(pl.col(end_column).cast(pl.List(pl.Int64)))
            proteins_df = proteins_df.with_columns(pl.col(start_column).cast(pl.List(pl.Int64)))
            proteins_df = proteins_df.with_columns(pl.col('peptide_index').cast(pl.List(pl.List(pl.Int64))))
            proteins_df = proteins_df.with_columns(pl.col('sample').cast(pl.List(pl.String)))
            proteins_df = proteins_df.with_columns(pl.col('condition').cast(pl.List(pl.String)))
            proteins_df = proteins_df.with_columns(pl.struct(start_column, end_column, seq_column, protacc_column, 'peptide_index', sample_column, condition_column).map_batches(lambda x: pl.Series(group_repetitive(x.struct.field(start_column), x.struct.field(end_column), x.struct.field(seq_column), x.struct.field(protacc_column), x.struct.field('peptide_index'), x.struct.field(sample_column), x.struct.field(condition_column))), return_dtype=pl.String).str.split('|').alias('repetitive'))
            proteins_df = proteins_df.with_columns(pl.col('repetitive').list.get(0).alias('start'))
            proteins_df = proteins_df.with_columns(pl.col('repetitive').list.get(1).alias('end'))
            proteins_df = proteins_df.with_columns(pl.col('repetitive').list.get(2).alias('peptide_index'))
            proteins_df = proteins_df.with_columns(pl.col('repetitive').list.get(3).alias('sequence'))
            proteins_df = proteins_df.with_columns(pl.col('repetitive').list.get(4).alias('sample'))
            proteins_df = proteins_df.with_columns(pl.col('repetitive').list.get(5).alias('condition'))
            proteins_df = proteins_df.with_columns((pl.col(start_column).str.split(delimiter)).alias(start_column))
            proteins_df = proteins_df.with_columns((pl.col(end_column).str.split(delimiter)).alias(end_column))
            proteins_df = proteins_df.with_columns((pl.col('peptide_index').str.split(delimiter)).alias('peptide_index'))
            proteins_df = proteins_df.with_columns((pl.col('sequence').str.split(delimiter)).alias('sequence'))
            proteins_df = proteins_df.with_columns((pl.col('sample').str.split(delimiter)).alias('sample'))
            proteins_df = proteins_df.with_columns((pl.col('condition').str.split(delimiter)).alias('condition'))
            proteins_df = proteins_df.explode('start', 'end', 'peptide_index','sequence', 'sample', 'condition')
            proteins_df = proteins_df.group_by(protacc_column).agg(pl.col(seq_column), pl.col('start'), pl.col('end'), pl.col('peptide_index'), pl.col('sample'), pl.col('condition'))
            proteins_df = proteins_df.rename({protacc_column:'accession'})

    return proteins_df


def parse_input(evidence_file: str, seq_column: str, protacc_column: str, intensity_column: str, start_column: str, end_column: str, delimiter: str, proteome_dict: dict[str,str], mod_pattern: str, sample_column: str, condition_column: str) -> tuple[pd.DataFrame, int, float]:
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
        A tuple of a pandas dataframe, an integer and a float. The dataframe contains 
        one protein per row and all peptides mapped to that protein in the 
        peptides_df, with their start position, end position and intensity.
        The integer is the number of peptides that are removed due to their 
        accession not appearing in the proteome file. The float is the total intensity 
        of the evidence file. 
    """
    peptides_df = read_id_output(evidence_file, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, sample_column, condition_column)
    
    # get peptides/proteins with protein accessions that do not appear in the proteome
    peptides = pl.Series(peptides_df.with_columns((pl.col(protacc_column).list.filter(~pl.element().is_in(list(proteome_dict.keys())))).alias('removed')).select('removed')).to_list()
    n_removed_proteins = set(itertools.chain.from_iterable(peptides))

    # remove all peptides occurring multiple times in different modification and charge states
    peptides_df = peptides_df.with_columns((pl.col(seq_column).str.replace_all('\(.*?\)', '')).alias(seq_column))
    # remove peptides with protein accessions that do not appear in the proteome 
    if start_column and end_column:
        peptides_df = peptides_df.with_columns(pl.col(start_column).list.gather(pl.col(protacc_column).list.eval(pl.arg_where(~pl.element().is_in(n_removed_proteins)).alias(start_column))))
        peptides_df = peptides_df.with_columns(pl.col(end_column).list.gather(pl.col(protacc_column).list.eval(pl.arg_where(~pl.element().is_in(n_removed_proteins)).alias(end_column))))
        peptides_df = peptides_df.group_by([protacc_column, seq_column, sample_column, condition_column]).agg(pl.col(start_column).first(), pl.col(end_column).first(), pl.col('peptide_index'))
    else:
        peptides_df = peptides_df.group_by([protacc_column, seq_column, sample_column, condition_column]).agg(pl.col('peptide_index'))
    peptides_df = peptides_df.with_columns(pl.col(protacc_column).list.gather(pl.col(protacc_column).list.eval(pl.arg_where(~pl.element().is_in(n_removed_proteins)).alias(protacc_column))))
    # remove peptides that are not annotated with any proteome accession
    peptides_df = peptides_df.remove(pl.col(protacc_column).list.len() == 0)

    logger.info(f'Peptides mapped to the following {len(n_removed_proteins)} proteins were removed since the proteins do not appear in the proteome fasta file: {n_removed_proteins}.')
    protein_df = prot_pep_link(peptides_df, seq_column, protacc_column, intensity_column, start_column, end_column, proteome_dict, mod_pattern, delimiter, sample_column, condition_column)
    if intensity_column:
        total_intens = peptides_df[intensity_column].sum()
    else:
        total_intens = 0

    return protein_df, len(n_removed_proteins), total_intens