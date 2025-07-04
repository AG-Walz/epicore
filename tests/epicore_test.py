# import modules
import pandas as pd
from epicore_utils.modules.parse_input import parse_input, proteome_to_dict
from epicore_utils.modules.map_result import map_pep_core
from epicore_utils.modules.compute_cores import compute_consensus_epitopes

# define test files
evidence_path = 'tests/evidence_file.csv'
proteome_path = 'tests/proteome.fasta'
result_path = 'tests/result_file.csv'

# define params
seq_column = 'sequence'
protacc_column = 'accessions'
intensity_column = ''
start_column = ''
end_column = ''
delimiter = ';'
mod_pattern = ''
min_overlap = 7
max_step_size = 5
min_epi_length = 10

# run epicore on test files
proteome_dict = proteome_to_dict(proteome_path)
protein_df, n_removed_peps = parse_input(evidence_path, seq_column , protacc_column, intensity_column, start_column, end_column, delimiter, proteome_dict, mod_pattern)
protein_df = compute_consensus_epitopes(protein_df, min_overlap, max_step_size, min_epi_length, intensity_column, mod_pattern, proteome_dict)
pep_cores_mapping = map_pep_core(evidence_path,protein_df,seq_column,protacc_column,start_column,end_column,intensity_column,delimiter,mod_pattern, proteome_dict)

# sort result df
pep_cores_mapping = pep_cores_mapping.sort_values(by='sequence').reset_index(drop=True).astype(str)

# import and sort expected result
result_file = pd.read_csv('tests/result_file.csv')
result_file = result_file.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected result
def test_equal():
  assert pep_cores_mapping.equals(result_file)