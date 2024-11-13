import pandas as pd 
import re
import numpy as np 
from Bio import SeqIO


def read_id_output(id_output, output_type, seq_column, protacc_column, pep_pos=''):
    '''
    input:
        - id_output: file that contains identification output
        - output_type: either csv, tsv or xlsx
        - seq_column: header of column containing peptide sequence information 
        - protacc_column: header of column containing protein accession information 
        - delimiter: delimiter that separates entries in one column
        - pep_pos: optional parameter containing comma separated list of headers of columns containing start and end position of peptide in protein 
    output:
        - peptides_df: pandas dataframe containing the columns sequence and protein accession (and optional start and stop column) of input file
    '''
    if output_type == 'CSV':
        peptides_df = pd.read_csv(id_output, delimiter=',')
    elif output_type == 'TSV':
        peptides_df = pd.read_csv(id_output, delimiter='\t')
    elif output_type == 'XLSX':
        print('xlsx')
        peptides_df = pd.read_excel(id_output)

    if pep_pos != '':
        # read in start and end position if information is provided
        start = pep_pos.split(',')[0]
        end = pep_pos.split(',')[1]
        peptides_df = peptides_df[[protacc_column,seq_column, start, end]]     

    else:
        peptides_df = peptides_df[[protacc_column,seq_column]]
    
    peptides_df = peptides_df.astype(str)
    return peptides_df


def proteome_to_df(proteome):
    '''
    input:
        - proteome: fasta file with the proteome
    output:
        - pandas dataframe containing one protein accession and its corresponding sequence per row
    '''
    proteome_df = pd.DataFrame(columns=['accession', 'sequence'])
    proteome = SeqIO.parse(open(proteome),'fasta')
    for protein in proteome:
        protein_entry = {'accession':protein.id, 'sequence':str(protein.seq)}
        proteome_df.loc[len(proteome_df)] = protein_entry
    print(proteome_df)
    return proteome_df


def get_prot_seq(accession, proteome_df):
    '''
    input:
        - accession (str)
        - reference proteome (fasta)
    output:
        - protein sequence corresponding to accession 
    '''
    #protein_seq = proteome_df.loc[proteome_df['accession'] == accession,'sequence'].iloc[0]
    if (proteome_df['accession'] == accession).any():
        protein_seq = proteome_df.loc[proteome_df['accession'] == accession,'sequence'].iloc[0]
    elif proteome_df['accession'].str.contains(accession).any():
        protein_seq = proteome_df.loc[proteome_df['accession'].str.contains(accession),'sequence'].iloc[0]
    else:
        raise Exception('The protein with {} does not occur in the given proteome! Please use the proteome that was used for the identification of the peptides.'.format(accession))
    return protein_seq    


def compute_pep_pos(peptide, prot_accession, proteome_df):
    '''
    input:
        - peptide: sequence of the peptide
        - prot_accession: accession of the protein 
        - proteome: path to file of proteome used for identification
    output:
        - start position(s) of peptide in protein
        - end position(s) of peptide in protein 
    '''    
    
    prot_seq = get_prot_seq(prot_accession, proteome_df)
    
    # regular expression for all occurrences of peptide in protein (also overlapping ones) 
    peptide_search = '(?=' + peptide + ')' 
    
    # get all occurrences of the peptide in the protein 
    pos_list = [[match.start(),match.start()+len(peptide)-1] for match in re.finditer(peptide_search,prot_seq)]
    if len(pos_list) == 0:
        raise Exception('The peptide {} does not occur in the protein with accession {} in the proteome you specified, but your input file provides evidence for that! Please use the proteome that was used for the identification of the peptides.'.format(peptide, prot_accession))
    
    pos_list = np.array(pos_list)
    start = pos_list[:,0]
    end = pos_list[:,1]

    return start, end


def group_repetitive(starts, ends):
    '''
    input:  
        - starts: list that contains all start positions of the peptide in the protein
        - end: list that contains all end positions of the peptide in the protein
    output:
        - if the the start and end positions only differ by one return the smallest start and highest end value
          else: starts and ends that were the input
    '''
    start_diffs = np.diff(starts)
    end_diffs = np.diff(ends)
    if set(start_diffs) == {1} and set(end_diffs) == {1}:
        return [min(starts)], [max(ends)]
    else:
        return starts, ends


def prot_pep_link(peptides_df, seq_column, protacc_column, delimiter, proteome, mod_delimiter, pep_pos=''):
    '''
    input:
        - peptides_df: pandas Dataframe, that contains one peptide per row
        - seq_column: header of the column in the input file containing the peptide sequences
        - protacc_column: header of the column in the input file containing the protein accessions
        - delimiter: delimiter that separates multiple values in one cell
        - proteome: path to the proteome used for identification 
        - mod_delimiter: comma separated string with delimiters for peptide_modifications
        - pep_pos: comma separated string with the header of the column in the input file containing the start and end positions of the peptides in the protein 
    output:
        - pandas DataFrame
            > row for each protein accession in input:
                protein accession 
                peptides matched 
                start positions of peptide in protein
                end positions of peptide in protein
    '''

    if pep_pos == '':
        # if the start and end positions of the peptides is not defined in the input evidence file 
        # 
        # load the proteome into a pandas dataframe
        proteome_df = proteome_to_df(proteome)

        proteins = pd.DataFrame(columns=[protacc_column, seq_column],dtype=object)
        proteins = pd.DataFrame(columns=['accession', 'sequence'])

        # create a row for each protein, that contains all peptides matched to the protein in the input evidence file
        for _, peptide in peptides_df.iterrows():
            prot_accessions = peptide[protacc_column].split(delimiter)
            for i, prot_accession in enumerate(prot_accessions):
                if prot_accession in proteins['accession'].values:
                    proteins.loc[proteins['accession'] == prot_accession, 'sequence'] = proteins.loc[proteins['accession'] == prot_accession, 'sequence'].apply(lambda x: x + [peptide[seq_column]])
                else:
                    protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]]}
                    proteins.loc[len(proteins)] = protein_entry

        # add the columns start and end
        proteins['start'] = [[] for _ in range(len(proteins))]
        proteins['end'] = [[] for _ in range(len(proteins))]

        # compute the peptide positions for each peptide of each protein
        for p, protein in proteins.iterrows():
            peptides = protein['sequence']
            accession = protein['accession']
            starts = []
            ends = []
            updated_peps = []

            for n_p, peptide in enumerate(peptides):
                updated_peps.append(peptide)

                # remove modifications in the sequence for the matching step 
                pattern = re.escape(mod_delimiter.split(',')[0]) + r'.*?' + re.escape(mod_delimiter.split(',')[1])
                peptide = re.sub(pattern,"",peptide)
                pep_start, pep_end = compute_pep_pos(peptide, accession, proteome_df)

                if pep_start.size == 0:
                    raise Exception('The peptide {} does not occur in the protein with accession {} in the proteome you specified, but your input file provides evidence for that! Please use the proteome that was used for the identification of the peptides.'.format(peptide, prot_accession))
                if pep_start.size > 1:
                    pep_start, pep_end = group_repetitive(pep_start,pep_end)
                    #print('CAUTION! The peptide sequence occurs multiple times')

                if len(pep_start) > 1:
                    for _ in pep_start[:-1]:
                        updated_peps.append(peptide)
                    
                    #print('CAUTION! The peptide sequence occurs multiple times at different positions')
                    #else:
                    print('CAUTION! The peptide sequence {} is part of a repetitive region and will be used as evidence of the entire repetitive region.'.format(peptide))
                    print(accession)
                    # TODO: repetitive region: write start and end position to evidence file and add to visualization
                
                # collect all start and end positions of the peptide in the protein
                for start in pep_start:
                    starts.append(start)
                for end in pep_end:
                    ends.append(end)
            
            # collect all start and end positions in the protein
            proteins.at[p, 'start'] = starts
            proteins.at[p, 'end'] = ends
            proteins.at[p, 'sequence'] = updated_peps
    else:
        # if the start and end positions of the peptides is defined in the input evidence file

        # get columns that contain start and end position 
        pep_start = pep_pos.split(',')[0]
        pep_end = pep_pos.split(',')[1]
        proteins = pd.DataFrame([protacc_column,seq_column, pep_start, pep_end],dtype=object)
        proteins = pd.DataFrame(columns=['accession', 'sequence', 'start','end'])

        for _, peptide in peptides_df.iterrows():

            # get all proteins associated with the peptide
            prot_accessions = peptide[protacc_column].split(delimiter)
            for i, prot_accession in enumerate(prot_accessions):
                if prot_accession in proteins['accession'].values:
                    proteins.loc[proteins['accession'] == prot_accession, 'sequence'] = proteins.loc[proteins['accession'] == prot_accession, 'sequence'].apply(lambda x: x + [peptide[seq_column]])
                    proteins.loc[proteins['accession'] == prot_accession, 'start'] = proteins.loc[proteins['accession'] == prot_accession, 'start'].apply(lambda x: x + [int(pos) for pos in [peptide[pep_start].split(delimiter)[i]]])
                    proteins.loc[proteins['accession'] == prot_accession, 'end'] = proteins.loc[proteins['accession'] == prot_accession, 'end'].apply(lambda x: x + [int(pos) for pos in [peptide[pep_end].split(delimiter)[i]]])
                else:
                    protein_entry = {'accession':prot_accession, 'sequence':[peptide[seq_column]], 'start':[int(pos) for pos in [peptide[pep_start].split(delimiter)[i]]], 'end':[int(pos) for pos in [peptide[pep_end].split(delimiter)[i]]]}
                    proteins.loc[len(proteins)] = protein_entry

    return proteins

def parse_input(input, input_type, seq_column, protacc_column, delimiter, proteome, mod_delimiter, pep_pos=''):
    '''
    input:
        - input: output file of search engine
        - input_type: CSV, TSV or XLSX
        - seq_column: header of the column in the input file containing the peptide sequences
        - protacc_column: header of the column in the input file containing the protein accessions
        - delimiter: delimiter that separates multiple values in one cell
        - proteome: path to the proteome used for identification 
        - mod_delimiter: comma separated string with delimiters for peptide_modifications
        - pep_pos: comma separated string with the header of the column in the input file containing the start and end positions of the peptides in the protein 
    output:
        - pandas DataFrame: 
            > row for each protein accession in input:
                protein accession 
                peptides matched 
                start positions of peptide in protein
                end positions of peptide in protein
    '''
    if pep_pos != '':
        peptides_df = read_id_output(input, input_type, seq_column, protacc_column,pep_pos)
        protein_df = prot_pep_link(peptides_df, seq_column, protacc_column, delimiter, proteome, mod_delimiter, pep_pos)
    else:
        peptides_df = read_id_output(input, input_type, seq_column, protacc_column)
        protein_df = prot_pep_link(peptides_df, seq_column, protacc_column, delimiter, proteome, mod_delimiter)
    return protein_df