"""
Reads in the evidence file and reference profile, computes the peptides positions in the proteome and links a protein accession with all peptides associated with the protein. 
"""

import pandas as pd 
import re
import numpy as np 
from Bio import SeqIO
import os 


def read_id_output(id_output, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter):
    '''
    input:
        - id_output: evidence file
        - seq_column: header of the column containing peptide sequence information in the evidence file
        - protacc_column: header of the column containing protein accession information in the evidence file
        - intensity_column: header of the column containing intensity information in the evidence file
        - delimiter: delimiter that separates multiple entries in one column in the evidence file
    output:
        - peptides_df: pandas dataframe containing the columns sequence, protein accession and peptide intensity (and optional start and stop column) of input file
    '''
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

    # check that the mandatory headers are provided and if all provided column headers are part of the evidence file
    if seq_column not in peptides_df.columns:
        if seq_column:
            raise Exception('The header {} does not exist in the provided evidence file. Please provide the correct header.'.format(seq_column))
        else:
            raise Exception('The header for the column containing the peptide sequences in the evidence file is mandatory. Please provide the correct header name.')
    if protacc_column not in peptides_df.columns:
        if protacc_column:
            raise Exception('The header {} does not exist in the provided evidence file. Please provide the correct header.'.format(protacc_column))
        else:
            raise Exception('The header for the column containing the protein accessions in the evidence file is mandatory. Please provide the correct header name.')
    if intensity_column not in peptides_df.columns:
        if intensity_column:
            raise Exception('The header {} does not exist in the provided evidence file. Please provide the correct header.'.format(intensity_column))
    if start_column not in peptides_df.columns:
        if start_column:
            raise Exception('The header {} does not exist in the provided evidence file. Please provide the correct header.'.format(start_column))
    if end_column not in peptides_df.columns:
        if end_column:
            raise Exception('The header {} does not exist in the provided evidence file. Please provide the correct header.'.format(end_column))

    # read in evidence file (start and stop position only when provided)
    if start_column and end_column:
        if intensity_column:
            peptides_df = peptides_df[[protacc_column,seq_column, intensity_column, start_column, end_column]]  
        else:
            peptides_df = peptides_df[[protacc_column,seq_column, start_column, end_column]]
        peptides_df[start_column] = peptides_df[start_column].apply(lambda x: list(x.split(delimiter)))
        peptides_df[end_column] = peptides_df[end_column].apply(lambda x: list(x.split(delimiter)))
        peptides_df[protacc_column] = peptides_df[protacc_column].apply(lambda x: list(x.split(delimiter)))
    else:
        if intensity_column:
            peptides_df = peptides_df[[protacc_column,seq_column, intensity_column]]
        else:
            peptides_df = peptides_df[[protacc_column,seq_column]]
        peptides_df[protacc_column] = peptides_df[protacc_column].apply(lambda x: list(set(x.split(delimiter))))

    # TODO: check for side effects when removing this statement
    #peptides_df = peptides_df.astype(str)

    # split accessions if peptide occurs in more than one protein 
    #peptides_df[protacc_column] = peptides_df[protacc_column].apply(lambda x: list(x.split(delimiter)))
    return peptides_df


def proteome_to_df(proteome):
    '''
    input:
        - proteome: fasta file containing the proteome
    output:
        - proteome_df: pandas dataframe containing one protein accession and its corresponding sequence per row
    '''
    proteome_df = pd.DataFrame(columns=['accession', 'sequence'])
    proteome = SeqIO.parse(open(proteome),'fasta')
    for protein in proteome:
        protein_entry = {'accession':protein.id, 'sequence':str(protein.seq)}
        proteome_df.loc[len(proteome_df)] = protein_entry
    return proteome_df


def get_prot_seq(accession, proteome_df):
    '''
    input:
        - accession: accession of protein
        - proteome_df: reference proteome used for the identification of the peptide as a pandas dataframe
    output:
        - protein_seq: protein sequence corresponding to accession 
    '''
    if (proteome_df['accession'] == accession).any():
        protein_seq = proteome_df.loc[proteome_df['accession'] == accession,'sequence'].iloc[0]
    else:
        raise Exception('The protein with accession "{}" does not occur in the given proteome! Please use the proteome that was used for the identification of the peptides.'.format(accession))
    return protein_seq    


def compute_pep_pos(peptide, accession, proteome_df):
    '''
    input:
        - peptide: sequence of the peptide
        - accession: accession of the protein 
        - proteome_df: reference proteome used for the identification of the peptide as a pandas dataframe
    output:
        - start: start position(s) of the peptide in the protein
        - end: end position(s) of the peptide in the protein 
    '''    
    
    prot_seq = get_prot_seq(accession, proteome_df)

    # regular expression for all occurrences of peptide in protein (also overlapping ones) 
    peptide_search = '(?=' + peptide + ')' 
    
    # get all occurrences of the peptide in the protein 
    pos_list = [[match.start(),match.start()+len(peptide)-1] for match in re.finditer(peptide_search,prot_seq)]
    if len(pos_list) == 0:
        raise Exception('The peptide {} does not occur in the protein with accession "{}"! Please use the proteome that was used for the identification of the peptides.'.format(peptide, accession))
    
    pos_list = np.array(pos_list)
    start = pos_list[:,0]
    end = pos_list[:,1]

    return start, end

def group_repetitive(starts, ends, peptide, accession):
    '''
    input:  
        - starts: list of all start positions of a peptide in a protein
        - ends: list of all end positions of a peptide in a protein
        - peptide: peptide sequence
        - accession: protein accession
    output:
        - starts: input start positions, but only the smallest start position for each repetitive region
        - ends: input end positions, but only the smallest end position for each repetitive region
    '''
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
        # TODO: find a better solution to inform the user 
        print('CAUTION! The peptide sequence {} is part of repetitive region(s) in protein {} and will be used as evidence of the entire repetitive region.'.format(peptide, accession))        
        return updated_starts, updated_ends
    else:
        return starts, ends


def prot_pep_link(peptides_df, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, proteome, mod_delimiter):
    '''
    input:
        - peptides_df: pandas Dataframe, that contains one peptide per row
        - seq_column: header of the column containing peptide sequence information in the evidence file
        - protacc_column: header of the column containing protein accession information in the evidence file
        - intensity_column: header of the column containing intensity information in the evidence file
        - delimiter: delimiter that separates multiple entries in one column in the evidence file
        - proteome: reference proteome used for the identification of the peptide as a pandas dataframe
        - mod_delimiter: comma separated string with delimiters for peptide modifications
    output:
        - proteins: pandas DataFrame
            > row for each protein accession in peptides_df:
                protein accession 
                peptides matched to the protein
                start positions of peptide in protein
                end positions of peptide in protein
                measured intensity of the peptide 
    '''

    if not start_column and not end_column:
        # if the start and end positions of the peptides is not defined in the input evidence file 
        
        # load the proteome into a pandas dataframe
        proteome_df = proteome_to_df(proteome)

        if intensity_column:
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'intensity'])
        else:
            proteins = pd.DataFrame(columns=['accession', 'sequence'])

        # create a row for each protein, that contains all peptides matched to the protein in the input evidence file
        for _, peptide in peptides_df.iterrows():
            # get all proteins matched to the peptide
            prot_accessions = peptide[protacc_column]
            for i, prot_accession in enumerate(prot_accessions):
                # TODO: remove this line! only for testing 
                if 'DECOY' not in prot_accession:
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
                intensity = protein['intensity']

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
                    updated_intens.append(intensity[n_p])
                
                # remove modifications in the sequence for the matching step (modifications are separated from the peptide by (), [] or by the user defined mod_delimiter)
                pattern = r'\(.*?\)'
                peptide = re.sub(pattern,"",peptide)
                pattern = r'\[.*?\]'
                peptide = re.sub(pattern,"",peptide)
                pattern = re.escape(mod_delimiter.split(',')[0]) + r'.*?' + re.escape(mod_delimiter.split(',')[1])
                peptide = re.sub(pattern,"",peptide)
                pep_start, pep_end = compute_pep_pos(peptide, accession, proteome_df)

                if pep_start.size == 0:
                    raise Exception('The peptide {} does not occur in the protein with accession {} in the proteome you specified, but your input file provides evidence for that! Please use the proteome that was used for the identification of the peptides.'.format(peptide, prot_accession))
                if pep_start.size > 1:
                    pep_start, pep_end = group_repetitive(pep_start,pep_end, peptide, accession)
                   
                if len(pep_start) > 1:
                    # if a peptide occurs multiple times in a protein
                    for _ in pep_start[:-1]:
                        updated_peps.append(peptide)
                        if intensity_column:
                            updated_intens.append(intensity[n_p])
                    print('The peptide sequence {} occurs multiple times in {}. It will be used as evidence for all occurrences.'.format(peptide, accession))
                
                # collect all start and end positions of the peptide in the protein
                for start in pep_start:
                    starts.append(start)
                for end in pep_end:
                    ends.append(end)
            
            proteins.at[p, 'start'] = starts
            proteins.at[p, 'end'] = ends
            proteins.at[p, 'sequence'] = updated_peps
            if intensity_column:
                proteins.at[p,'intensity'] = updated_intens
    else:
        # if the start and end positions of the peptides is defined in the input evidence file

        if intensity_column:
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'intensity', 'start','end'])
        else:
            proteins = pd.DataFrame(columns=['accession', 'sequence', 'start','end'])

        for _, peptide in peptides_df.iterrows():

            # get all proteins associated with the peptide
            prot_accessions = peptide[protacc_column]
            for i, prot_accession in enumerate(prot_accessions):

                if prot_accession in proteins['accession'].values:
                    # create new row for accessions seen before
                    proteins.loc[proteins['accession'] == prot_accession, 'sequence'] = proteins.loc[proteins['accession'] == prot_accession, 'sequence'].apply(lambda x: x + [peptide[seq_column]])
                    if intensity_column:
                        proteins.loc[proteins['accession'] == prot_accession, 'intensity'] = proteins.loc[proteins['accession'] == prot_accession, 'intensity'].apply(lambda x: x + [peptide[intensity_column]])
                    proteins.loc[proteins['accession'] == prot_accession, 'start'] = proteins.loc[proteins['accession'] == prot_accession, 'start'].apply(lambda x: x + [peptide[start_column][i]])
                    proteins.loc[proteins['accession'] == prot_accession, 'end'] = proteins.loc[proteins['accession'] == prot_accession, 'end'].apply(lambda x: x + [peptide[end_column][i]])
                    
                else:
                    # create new row for accessions not seen before
                    if intensity_column:
                        protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]], 'intensity':[peptide[intensity_column]], 'start':[int(pos) for pos in [peptide[start_column][i]]], 'end':[int(pos) for pos in [peptide[end_column][i]]]}
                    else:
                        protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]], 'start':[peptide[start_column][i]], 'end':[peptide[end_column][i]]}
                    
                    proteins.loc[len(proteins)] = protein_entry

    return proteins

def parse_input(input, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, proteome, mod_delimiter):
    '''
    input:input:
        - peptides_df: pandas Dataframe, that contains one peptide per row
        - seq_column: header of the column containing peptide sequence information in the evidence file
        - protacc_column: header of the column containing protein accession information in the evidence file
        - intensity_column: header of the column containing intensity information in the evidence file
        - delimiter: delimiter that separates multiple entries in one column in the evidence file
        - proteome: reference proteome used for the identification of the peptide as a pandas dataframe
        - mod_delimiter: comma separated string with delimiters for peptide modifications
    output:
        - pandas DataFrame: protein_df
            > row for each protein accession in peptides_df:
                protein accession 
                peptides matched to the protein
                start positions of peptide in protein
                end positions of peptide in protein
                measured intensity of the peptide 
    '''
    
    peptides_df = read_id_output(input, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter)
    protein_df = prot_pep_link(peptides_df, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, proteome, mod_delimiter)
    return protein_df