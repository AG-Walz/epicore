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

def run_epicore(proteome_path:str, evidence_path: str, seq_column: str, protacc_column: str, intensity_column: str, start_column: str, end_column: str, delimiter: str, mod_pattern: str, min_overlap: int, max_step_size: int, min_epi_length: int) -> pd.DataFrame:
  """Run epicore on an evidence file.

  Args:
    proteome_path: Path to the proteome fasta.
    evidence_path: Path to the evidence file.
    seq_column: The header of the column containing peptide sequence
        information in the evidence file.
    protacc_column: The header of the column containing protein accession
        information in the evidence file.
    intensity_column: The header of the column containing intensity
        information in the evidence file.
    start_column: The header of the column containing the start positions
        of the peptides in the protein sequences.
    end_column: The header of the column containing the end position of
        the peptides in the protein sequences.
    delimiter: The delimiter that separates multiple entries in one column 
        in the evidence file.
    mod_pattern: A comma separated string with delimiters for peptide
        modifications.
    min_overlap: An integer of the minimal overlap between two epitopes
        to be grouped to the same consensus epitope.
    max_step_size: An integer of the maximal distance between the start
        position of two epitopes to be grouped to the same consensus
        epitope.
    min_epi_len: An integer of the minimal length of a consensus epitope.

  Returns:
      The sorted evidence_file with four additional columns containing the
      whole and core sequence and total and relative intensity of each
      consensus epitope group, to which the peptide of the row belongs.
  """
  proteome_dict = proteome_to_dict(proteome_path)
  protein_df, _, total_intens = parse_input(evidence_path, seq_column , protacc_column, intensity_column, start_column, end_column, delimiter, proteome_dict, mod_pattern)
  protein_df = compute_consensus_epitopes(protein_df, min_overlap, max_step_size, min_epi_length, intensity_column, mod_pattern, proteome_dict, total_intens)
  pep_cores_mapping = map_pep_core(evidence_path,protein_df,seq_column,protacc_column,start_column,end_column,intensity_column,delimiter,mod_pattern, proteome_dict)
  pep_cores_mapping = pep_cores_mapping.sort_values(by='sequence').reset_index(drop=True).astype(str)
  return pep_cores_mapping


#############################
# Test one
#############################

# run epicore on test file one 
pep_cores_mapping_one = run_epicore(proteome_path, evidence_path, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length)

# read in expected result and sort it
result_file = pd.read_csv(path_resultone)
result_file = result_file.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_one():
  assert pep_cores_mapping_one.equals(result_file)




#############################
# Test two
#############################
max_step_size = 3
min_overlap = 8

# run epicore on test file two
pep_cores_mapping_two = run_epicore(proteome_path, evidence_path, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length)

# read in expected result and sort it
result_file_two = pd.read_csv(path_resulttwo)
result_file_two = result_file_two.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected result
def test_two():
  assert pep_cores_mapping_two.equals(result_file_two)