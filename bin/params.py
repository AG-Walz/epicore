'''
this file defines all parameters of localplateau that the user can change
'''

# put the path of your input file here
mhcquant_out = 'test_data/DN02_Brain_class2.tsv'

# put the path of the proteome used for identification of the peptides here
fasta_proteome = 'test_data/mod_proteome.fasta'

# put the minimal overlap of two peptides to still be grouped together here
min_overlap = 11

# put the maximal step size between two peptides to still be grouped together here
max_step_size = 5

# put the minimal length of the core epitopes here
min_epi_length = 9

# put the header of the column that contains the peptide sequences here
seq_column = 'sequence'

# put the header of the column that contains the protein accessions of the protein associated with the peptide in the identification step 
protacc_column = 'accessions'

# put the header of the column that contains the start position of the peptide in the protein (if that column does not exist use '')
start_column = ''

# put the header of the column that contains the end position of the peptide in the protein (if that column does not exist use '')
end_column = ''

# put the directory where you want to store your results here
out_dir =  'Brain_class2_mhcquant'

# comma separated string that contains the characters beginning and ending modifications in the peptide sequence 
mod_delimiter = '(,)' 

# put the delimiter that separates multiple entries in one column in the input file 
delimiter = ';'

# put a comma separated string of protein accessions for which you want to visualize the plateau result 
prot_accession = None

# put the path of a precomputed localplateau result here if you dont want to recompute everything and only visualize some proteins
plateau_csv = None