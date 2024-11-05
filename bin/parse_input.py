import pandas as pd 
from bin.compute_cores import get_prot_seg
import re
import numpy as np 

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
        
    return peptides_df

        
def compute_pep_pos(peptide, prot_accession, proteome):
    '''
    input:
        - peptide: sequence of the peptide
        - prot_accession: accession of the protein 
        - proteome: path to file of proteome used for identification
    output:
        - start position(s) of peptide in protein
        - end position(s) of peptide in protein 
    '''      
    prot_seq = get_prot_seg(prot_accession, proteome)
    pos_list = prot_seq.split(peptide)
    if len(pos_list) > 2:
        print('CAUTION! The peptide occurs multiple times.')
    start = []
    end = []
    current_pos = 0
    for pos in pos_list[:-1]:
        current_pos += len(pos)
        start.append(current_pos)
        end.append(current_pos+len(peptide)-1)
        current_pos += len(peptide)
    return start, end



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
        proteins = pd.DataFrame(columns=[protacc_column, seq_column],dtype=object)
        for _, peptide in peptides_df.iterrows():
            prot_accessions = peptide[protacc_column].split(delimiter)
            for i, prot_accession in enumerate(prot_accessions):
                if prot_accession in proteins[protacc_column].values:
                    proteins.loc[proteins[protacc_column] == prot_accession, seq_column] = proteins.loc[proteins[protacc_column] == prot_accession, seq_column].apply(lambda x: x + [peptide[seq_column]])
                else:
                    protein_entry = {protacc_column:prot_accession, seq_column:[peptide[seq_column]]}
                    proteins.loc[len(proteins)] = protein_entry
        proteins['start'] = [[] for _ in range(len(proteins))]
        proteins['end'] = [[] for _ in range(len(proteins))]
        for p, protein in proteins.iterrows():
            peptides = protein[seq_column]
            accession = protein[protacc_column]
            starts = []
            ends = []
            for _, peptide in enumerate(peptides):
                pattern = re.escape(mod_delimiter.split(',')[0]) + r'.*?' + re.escape(mod_delimiter.split(',')[1])
                peptide = re.sub(pattern,"",peptide)
                pep_start, pep_end = compute_pep_pos(peptide, accession, proteome)
                if pep_start == []:
                    print('CAUTION! The peptide sequence does not occur in the protein sequence.')
                    print(peptide)
                    print(accession)
                for start in pep_start:
                    starts.append(start)
                for end in pep_end:
                    ends.append(end)
            proteins.at[p, 'start'] = starts
            proteins.at[p, 'end'] = ends
    else:
        peptides_df = peptides_df[['accessions','sequence', 'start', 'end']]
        proteins = pd.DataFrame(columns=['accession', 'sequence', 'start','end'])

        for _, peptide in peptides_df.iterrows():
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