"""
Assigns each peptide in the evidence files its core epitopes, the total intensity of that core epitope and the relative core intensity. 
"""

import pandas as pd
import re
import os 

def read_entire_id_output(id_output):
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

def map_pep_core(evidence_file, protein_df, seq_column, protacc_column, start_column, end_column, intensity_column, delimiter, mod_delimiter):
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
        mod_delimiter: A comma separated string with delimiters for peptide
            modifications

    Returns:
        The evidence_file with four additional columns containing the whole and 
        core sequence and total and relative intensity of each consensus 
        epitope group, to which the peptide of the row belongs.

    Raises:
        Exception: If the mappings are contradictory.
    """
    # add the columns whole and core epitopes to the input evidence
    evidence_file_df = read_entire_id_output(evidence_file)
    evidence_file_df['whole_epitopes'] = [[] for _ in range(len(evidence_file_df))]
    evidence_file_df['core_epitopes'] = [[] for _ in range(len(evidence_file_df))]
    evidence_file_df['proteome_occurence'] = [[] for _ in range(len(evidence_file_df))]
    if intensity_column:
        evidence_file_df['total_core_intensity'] = [[] for _ in range(len(evidence_file_df))]
        evidence_file_df['relative_core_intensity'] = [[] for _ in range(len(evidence_file_df))]

    for r, row in evidence_file_df.iterrows():

        # loop over all proteins mapped to the peptide in the evidence file 
        for mapping, accession in enumerate(row[protacc_column].split(delimiter)):
            #TODO: only for testing!
            if 'DECOY' not in accession:
                
                sequence = row[seq_column]
                    
                # get protein row that contains the current peptide sequence and is associated with the protein from the evidence file
                prot_row = protein_df[(protein_df['accession'] == accession) & protein_df['sequence'].map(lambda x: sequence in x)]
                
                # indices of peptides that match the sequence of the peptide and the accession of the mapped protein
                idx = [i for i, x in enumerate(prot_row['sequence'].to_list()[0]) if x == sequence]
                
                if len(idx) > 1:
                    # check if multiple occurrence due to modification                        
                    wo_mod = [re.sub(r"[\[\(].*?[\]\)]","",prot_row['sequence'].to_list()[0][i]) for i in idx]
                    pattern = re.escape(mod_delimiter.split(',')[0]) + r'.*?' + re.escape(mod_delimiter.split(',')[1])
                    wo_mod = wo_mod + [re.sub(pattern,"",prot_row['sequence'].to_list()[0][i]) for i in idx]
                    if len(set(wo_mod)) > 1:
                        raise Exception('Please check your evidence file. There are peptides with different sequences mapped to the same position in the protein.')
                
                
                # get core and whole epitope associated with the peptide in the evidence file
                for id in idx:
                    mapped_group = prot_row['sequence_group_mapping'].to_list()[0][id]
                    evidence_file_df.at[r,'core_epitopes'].append(prot_row['consensus_epitopes'].to_list()[0][mapped_group])
                    evidence_file_df.at[r,'whole_epitopes'].append(prot_row['whole_epitopes'].to_list()[0][mapped_group])
                    if intensity_column:
                        evidence_file_df.at[r,'total_core_intensity'].append(str(prot_row['core_epitopes_intensity'].to_list()[0][mapped_group]))
                        evidence_file_df.at[r,'relative_core_intensity'].append(str(prot_row['relative_core_intensity'].to_list()[0][mapped_group]))
                    prot_occurence = accession +':'+ str(prot_row['core_epitopes_start'].to_list()[0][mapped_group]) + '-' + str(prot_row['core_epitopes_end'].to_list()[0][mapped_group])
                    evidence_file_df.at[r,'proteome_occurence'].append(prot_occurence)
            

    
        # convert list to delimiter separated strings
        evidence_file_df.at[r,'core_epitopes'] = delimiter.join(evidence_file_df.at[r,'core_epitopes'])
        evidence_file_df.at[r,'whole_epitopes'] = delimiter.join(evidence_file_df.at[r,'whole_epitopes'])
        if intensity_column:
            evidence_file_df.at[r,'total_core_intensity'] = delimiter.join(evidence_file_df.at[r,'total_core_intensity'])
            evidence_file_df.at[r,'relative_core_intensity'] = delimiter.join(evidence_file_df.at[r,'relative_core_intensity'])
        evidence_file_df.at[r,'proteome_occurence'] = delimiter.join(evidence_file_df.at[r,'proteome_occurence'])
        
    return evidence_file_df


