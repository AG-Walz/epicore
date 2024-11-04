import pandas as pd
import re

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
                print(accession)
                print([prot_row['sequence'].to_list()[0][i] for i in idx])
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