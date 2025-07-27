# import modules
import pandas as pd
from epicore_utils.modules.parse_input import parse_input, proteome_to_dict
from epicore_utils.modules.map_result import map_pep_core
from epicore_utils.modules.compute_cores import compute_consensus_epitopes

# define test files
evidence_path = 'tests/evidence_file.csv'
proteome_path = 'tests/proteome.fasta'
path_resultone = 'tests/result_one.csv'
path_resulttwo = 'tests/result_two.csv'

# define params
seq_column = 'sequence'
protacc_column = 'accessions'
intensity_column = 'intensity'
start_column = ''
end_column = ''
delimiter = ';'
mod_pattern = ''
min_overlap = 7
max_step_size = 5
min_epi_length = 10

# run epicore on test files
proteome_dict = proteome_to_dict(proteome_path)
protein_df, n_removed_peps, total_intens = parse_input(evidence_path, seq_column , protacc_column, intensity_column, start_column, end_column, delimiter, proteome_dict, mod_pattern)
protein_df = compute_consensus_epitopes(protein_df, min_overlap, max_step_size, min_epi_length, intensity_column, mod_pattern, proteome_dict, total_intens)
pep_cores_mapping = map_pep_core(evidence_path,protein_df,seq_column,protacc_column,start_column,end_column,intensity_column,delimiter,mod_pattern, proteome_dict)

# sort result df
pep_cores_mapping = pep_cores_mapping.sort_values(by='sequence').reset_index(drop=True).astype(str)
result_file = pd.read_csv(path_resultone)
result_file = result_file.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_one():
  assert pep_cores_mapping.equals(result_file)

max_step_size = 3
min_overlap = 8

# run epicore on test files
protein_df, n_removed_peps, total_intens = parse_input(evidence_path, seq_column , protacc_column, intensity_column, start_column, end_column, delimiter, proteome_dict, mod_pattern)
protein_df = compute_consensus_epitopes(protein_df, min_overlap, max_step_size, min_epi_length, intensity_column, mod_pattern, proteome_dict, total_intens)
pep_cores_mapping_two = map_pep_core(evidence_path,protein_df,seq_column,protacc_column,start_column,end_column,intensity_column,delimiter,mod_pattern, proteome_dict)

# sort result df
pep_cores_mapping_two = pep_cores_mapping_two.sort_values(by='sequence').reset_index(drop=True).astype(str)
result_file_large = pd.read_csv(path_resulttwo)
result_file_large = result_file_large.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected result
def test_two():
  assert pep_cores_mapping_two.equals(result_file_large)