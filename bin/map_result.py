import pandas as pd
import re

def map_pep_core(input_tsv, protein_df, delimiter, position_boolean=False):
    '''
    input:
        - input_tsv: MHCquant output
        - protein_df: pandas DataFrame with one protein per row and all core and whole epitopes matched to that position 
    output:
        - input_tsv with two more columns, one for the whole epitope and one for the core epitope as a pandas dataframe
    '''
    
    # add the columns whole and core epitopes to the input evidence
    MHCquant_out = pd.read_csv(input_tsv, delimiter='\t')
    MHCquant_out['whole_epitopes'] = [[] for _ in range(len(MHCquant_out))]
    MHCquant_out['core_epitopes'] = [[] for _ in range(len(MHCquant_out))]
    
    for r, row in MHCquant_out.iterrows():

        seen_comps = {}
        idx_pos = 0

        # loop over all proteins mapped to the peptide in the evidence file 
        for mapping, accession in enumerate(row['accessions'].split(';')):
            
            # reset index position for each new mapping 
            idx_pos = 0

            if position_boolean:
                # get start and end index of that peptide in the protein
                start = row['start'].split(delimiter)[mapping]
                end = row['end'].split(delimiter)[mapping]
                sequence = row['sequence']

                # get protein row that contains the current peptide sequence and is associated with the protein from the evidence file
                prot_row = protein_df[(protein_df['accession'] == accession) & protein_df['sequence'].map(lambda x: sequence in x)]


                # indices of peptides that match the sequence of the peptide, the accession of the mapped protein and the start and end position in the protein
                idx = [i for i, x in enumerate(zip(prot_row['start'].to_list()[0],prot_row['end'].to_list()[0])) if (x[0] == int(start) and x[1] == int(end))]

                #if str(seq)

                if len(idx) > 1: 
                    # check if multiple occurrence due to modification
                    wo_mod = [re.sub(r"\(.*?\)","",prot_row['sequence'].to_list()[0][i]) for i in idx]
                    if len(set(wo_mod)) > 1:
                        raise Exception('Please check your evidence file. There are peptides with different sequences mapped to the same position in the protein.')
                
                
            else:
                
                # TODO check if all evidence used
                sequence = row['sequence']


                if str(sequence) + str(accession) not in seen_comps:
                    seen_comps[str(sequence) + str(accession)] = 1
                else:
                    idx_pos = seen_comps[str(sequence) + str(accession)]
                    seen_comps[str(sequence) + str(accession)] += 1


                # get protein row that contains the current peptide sequence and is associated with the protein from the evidence file
                prot_row = protein_df[(protein_df['accession'] == accession) & protein_df['sequence'].map(lambda x: sequence in x)]
                
                # indices of peptides that match the sequence of the peptide and the accession of the mapped protein
                idx = [i for i, x in enumerate(prot_row['sequence'].to_list()[0]) if x == sequence]

                if len(idx) > 1: 
                    # check if multiple occurrence due to modification
                    wo_mod = [re.sub(r"\(.*?\)","",prot_row['sequence'].to_list()[0][i]) for i in idx]
                    if len(set(wo_mod)) > 1:
                        raise Exception('Please check your evidence file. There are peptides with different sequences mapped to the same position in the protein.')
                    
                '''if len(idx) > 1:
                    print(sequence)
                    print(accession)
                    for i in idx:
                        mapped_group = prot_row['sequence_group_mapping'].to_list()[0][i]
                        print(mapped_group)
                        print(prot_row['consensus_epitopes'].to_list()[0][mapped_group]) 
                    print('__________________')'''
                
            # get core and whole epitope associated with the peptide in the evidence file
            if idx_pos > 0:
                print(sequence)
                print(accession)
                print(idx_pos)
                print(prot_row['sequence_group_mapping'].to_list()[0][idx[idx_pos]])
                print(seen_comps)
                print('_____________')

            mapped_group = prot_row['sequence_group_mapping'].to_list()[0][idx[idx_pos]]

            MHCquant_out.at[r,'core_epitopes'].append(prot_row['consensus_epitopes'].to_list()[0][mapped_group])
            MHCquant_out.at[r,'whole_epitopes'].append(prot_row['whole_epitopes'].to_list()[0][mapped_group])

    
        # convert list to ; separated strings
        MHCquant_out.at[r,'core_epitopes'] = ';'.join(MHCquant_out.at[r,'core_epitopes'])
        MHCquant_out.at[r,'whole_epitopes'] = ';'.join(MHCquant_out.at[r,'whole_epitopes'])
        
    return MHCquant_out


