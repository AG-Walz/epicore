"""
Assigns each peptide in the evidence files its core epitopes, the total intensity of that core epitope and the relative core intensity. 
"""

import pandas as pd
import re
import os 
from itertools import repeat
import logging
logger = logging.getLogger(__name__)

def read_entire_id_output(id_output: str) -> pd.DataFrame:
    """Read in the entire evidence file.
    
    Args:
        id_output: The string of the path to the evidence file.
    
    Returns:
        A pandas dataframe containing the evidence file.

    Raises:
        Exception: If the file type of the provided evidence file is not 
            supported.
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
        raise Exception('The file type of your evidence file is not supported. Please use an evidence file that has one of the following file types: csv, tsv, xlsx')
    return peptides_df

def map_pep_core(evidence_file: str, protein_df: pd.DataFrame, seq_column: str, protacc_column: str, start_column: str, end_column: str, intensity_column: str, delimiter: str, mod_pattern: str, proteome_dict: dict[str,str]) -> pd.DataFrame:
    """Map computed consensus epitope groups to the input evidence_file.
    
    Args:
        evidence_file: The string of the path to the evidence file.
        protein_df: A pandas dataframe containing one protein per row.
        seq_column: The string of the header of the column containing 
            peptide sequence information in the evidence file.
        protacc_column: The string of the header of the column containing 
            protein accession information in the evidence file.
        start_column: The string of the header of the column containing the 
            start positions of peptides in proteins.
        end_column: The string of the header of the column containing the end 
            position of peptides in proteins.
        intensity_column: The string of the header of the column containing 
            intensity information in the evidence file.
        delimiter: The delimiter that separates multiple entries in one column 
            in the evidence file.
        mod_pattern: A comma separated string with delimiters for peptide
            modifications

    Returns:
        The evidence_file with four additional columns containing the whole and 
        core sequence and total and relative intensity of each consensus 
        epitope group, to which the peptide of the row belongs.

    Raises:
        Exception: If the mappings are contradictory.
    """

    # remove accessions from evidence file that do not occur in the proteome 
    evidence_file_df = read_entire_id_output(evidence_file)
    idx_cols = evidence_file_df.columns[evidence_file_df.map(lambda cell: delimiter in str(cell)).any()]
    evidence_file_df[idx_cols] = evidence_file_df[idx_cols].map(lambda cell: cell.split(delimiter))
    evidence_file_df['indices'] = evidence_file_df.apply(lambda row: [idx for idx, prot in enumerate(row[protacc_column]) if prot in proteome_dict.keys()], axis=1)
    evidence_file_df.apply(lambda row: [[row[col][idx] for idx in row['indices']] for col in idx_cols], axis=1).to_csv('evidence_file.csv')
    evidence_file_df[idx_cols] = pd.DataFrame(evidence_file_df.apply(lambda row: [[row[col][idx] for idx in row['indices']] for col in idx_cols], axis=1).to_list(), index=evidence_file_df.index)

    # add the columns whole and core epitopes to the input evidence
    evidence_file_df = evidence_file_df.assign(
        whole_epitopes = list(repeat([], len(evidence_file_df))),
        core_epitopes = list(repeat([], len(evidence_file_df))),
        proteome_occurrence = list(repeat([], len(evidence_file_df)))
    )
    if intensity_column:
        evidence_file_df = evidence_file_df.assign(
            total_core_intensity = list(repeat([], len(evidence_file_df))),
            relative_core_intensity = list(repeat([], len(evidence_file_df)))
        )

    for r, row in evidence_file_df.iterrows():

        # loop over all proteins mapped to the peptide in the evidence file 
        for mapping, accession in enumerate(row[protacc_column]):
            sequence = row[seq_column]
                    
            # get protein row that contains the current peptide sequence and is associated with the protein from the evidence file
            prot_row = protein_df[(protein_df['accession'] == accession) & protein_df['sequence'].map(lambda x: sequence in x)]
                
            # indices of peptides that match the sequence of the peptide and the accession of the mapped protein
            idx = [i for i, seq in enumerate(prot_row['sequence'].iloc[0]) if seq == sequence]
                
            if len(idx) > 1:
                # check if multiple occurrence due to modification                        
                wo_mod = [re.sub(r"[\[\(].*?[\]\)]","",prot_row['sequence'].iloc[0][i]) for i in idx]
                if mod_pattern:
                    pattern = re.escape(mod_pattern.split(',')[0]) + r'.*?' + re.escape(mod_pattern.split(',')[1])
                    wo_mod = wo_mod + [re.sub(pattern,"",prot_row['sequence'].iloc[0][i]) for i in idx]
                if len(set(wo_mod)) > 1:
                    raise Exception('Please check your evidence file. There are peptides with different sequences mapped to the same position in the protein.')
                
                
            # get core and whole epitope associated with the peptide in the evidence file
            for id in idx:
                mapped_group = prot_row['sequence_group_mapping'].iloc[0][id]
                evidence_file_df.at[r,'core_epitopes'].append(prot_row['consensus_epitopes'].iloc[0][mapped_group])
                evidence_file_df.at[r,'whole_epitopes'].append(prot_row['whole_epitopes'].iloc[0][mapped_group])
                if intensity_column:
                    evidence_file_df.at[r,'total_core_intensity'].append(str(prot_row['core_epitopes_intensity'].iloc[0][mapped_group]))
                    evidence_file_df.at[r,'relative_core_intensity'].append(str(prot_row['relative_core_intensity'].iloc[0][mapped_group]))
                prot_occurence = accession +':'+ str(prot_row['core_epitopes_start'].iloc[0][mapped_group]) + '-' + str(prot_row['core_epitopes_end'].iloc[0][mapped_group])
                evidence_file_df.at[r,'proteome_occurrence'].append(prot_occurence)

        
        # convert list to delimiter separated strings
        evidence_file_df.at[r,'core_epitopes'] = delimiter.join(evidence_file_df.at[r,'core_epitopes'])
        evidence_file_df.at[r,'whole_epitopes'] = delimiter.join(evidence_file_df.at[r,'whole_epitopes'])
        if intensity_column:
            evidence_file_df.at[r,'total_core_intensity'] = delimiter.join(evidence_file_df.at[r,'total_core_intensity'])
            evidence_file_df.at[r,'relative_core_intensity'] = delimiter.join(evidence_file_df.at[r,'relative_core_intensity'])
        evidence_file_df.at[r,'proteome_occurrence'] = delimiter.join(evidence_file_df.at[r,'proteome_occurrence'])

    evidence_file_df = evidence_file_df.drop(columns='indices') 
    return evidence_file_df


def gen_epitope_df(protein_df: pd.DataFrame) -> pd.DataFrame:
    """Generate dataframe that has one epitope per row.
    
    Args:
        protein_df: A pandas dataframe containing one protein per row.

    Returns:
        A reordered version of protein_df were each row stores one epitope.
    """
    # include intensity columns if present
    if ('core_epitopes_intensity' not in protein_df.columns) and ('relative_core_intensity' not in protein_df.columns):
        cols = ['whole_epitopes', 'consensus_epitopes','landscape', 'grouped_peptides_sequence']
    else:
        cols = ['whole_epitopes', 'consensus_epitopes','landscape', 'grouped_peptides_sequence', 'relative_core_intensity', 'core_epitopes_intensity']

    cols_acc = cols + ['accession']

    # separate each epitope in one row
    protein_df_long = protein_df.explode(cols)
    protein_df_long = protein_df_long.astype(str)
    epitopes_grouped_df = protein_df_long[cols_acc].groupby(cols)
    epitopes_grouped_df = epitopes_grouped_df.agg({'accession':lambda x:','.join(x)}).reset_index()
    
    logger.info(f'{len(epitopes_grouped_df)} unique epitopes were computed.')

    return epitopes_grouped_df
