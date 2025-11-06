# import modules
import pandas as pd
from epicore_utils.modules.parse_input import parse_input, proteome_to_dict
from epicore_utils.modules.map_result import map_pep_core
from epicore_utils.modules.compute_cores import compute_consensus_epitopes
import polars as pl



def run_epicore(proteome_path:str, evidence_path: str, seq_column: str, protacc_column: str, intensity_column: str, start_column: str, end_column: str, delimiter: str, mod_pattern: str, min_overlap: int, max_step_size: int, min_epi_length: int, sample_column) -> pd.DataFrame:
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
  protein_df, _, total_intens = parse_input(evidence_path, seq_column , protacc_column, intensity_column, start_column, end_column, delimiter, proteome_dict, mod_pattern, sample_column)
  protein_df = protein_df.to_pandas()
  protein_df = compute_consensus_epitopes(protein_df, min_overlap, max_step_size, min_epi_length, intensity_column, mod_pattern, proteome_dict, total_intens)
  pep_cores_mapping = map_pep_core(evidence_path,protein_df,seq_column,protacc_column,start_column,end_column,intensity_column,delimiter,mod_pattern, proteome_dict)
  pep_cores_mapping = pep_cores_mapping.sort_values(by='sequence').reset_index(drop=True).astype(str)
  return pep_cores_mapping

# define test files
evidence_path = 'tests/evidence_file.csv'
#evidence_path_two = 'tests/evidence_file.csv'
evidence_path_three = 'tests/evidence_file_cohort.csv'
evidence_path_four = 'tests/evidence_file_modifications.csv'
evidence_path_five = 'tests/evidence_file_five.csv'
large_evidence = 'tests/UPN52.tsv'
path_result_one = 'tests/result_one.csv'
path_result_two = 'tests/result_two.csv'
path_result_three = 'tests/result_three.csv'
path_result_four = 'tests/results_four.csv'
path_result_five = 'tests/result_five.csv'
large_result = 'tests/UPN52_cw.tsv'

large_fasta = 'tests/spHUMANwoi_130927_CLL_mut.fasta'
proteome_path = 'tests/proteome.fasta'

# define params
seq_column = 'sequence'
protacc_column = 'accessions'
intensity_column = None
sample_column = 'sample'
start_column = 'start'
end_column = 'end'
delimiter = ';'
mod_pattern = ''
min_overlap = 7
max_step_size = 4
min_epi_length = 10


#############################
# Test one (one sample)
#############################

# run epicore on test file one 
pep_cores_mapping_one = run_epicore(proteome_path, evidence_path, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column)

# read in expected result and sort it
result_file_one = pd.read_csv(path_result_one)
result_file_one = result_file_one.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_one():
  assert pep_cores_mapping_one.equals(result_file_one)
  
  
#############################
# Test two (multiple samples)
#############################

pep_cores_mapping_two = run_epicore(proteome_path, evidence_path_three, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column)

# read in expected result and sort it
result_file_two = pd.read_csv(path_result_three)
result_file_two = result_file_two.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_two():
  assert pep_cores_mapping_two.equals(result_file_two)


#############################
# Test three (multiple samples, multiple charges, modifications)
#############################

min_epi_length = 3

# run epicore on test file one 
pep_cores_mapping_three = run_epicore(proteome_path, evidence_path_four, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column)

# read in expected result and sort it
result_file_three = pd.read_csv(path_result_four)
result_file_three = result_file_three.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_three():
  assert pep_cores_mapping_three.equals(result_file_three)


#############################
# Test four 
#############################

max_step_size = 5
min_epi_length = 3

# run epicore on test file one 
pep_cores_mapping_four = run_epicore(proteome_path, evidence_path_five, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column)
pep_cores_mapping_four = pep_cores_mapping_four.sort_values(by='sequence')

# read in expected result and sort it
result_file_four = pd.read_csv(path_result_five)
result_file_four = result_file_four.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_four():
  assert pep_cores_mapping_four.equals(result_file_four)

#############################
# Test five (one sample, large test)
#############################

max_step_size = 5
min_epi_length = 13
min_overlap = 11

# run epicore on test file one 
pep_cores_mapping_five = run_epicore(large_fasta, large_evidence, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column)
pep_cores_mapping_five = pep_cores_mapping_five.sort_values(by='sequence')
pep_cores_mapping_five = pep_cores_mapping_five.loc[:, ~pep_cores_mapping_five.columns.str.contains('^Unnamed')]
pep_cores_mapping_five['start'] = pep_cores_mapping_five['start'].str.split(';')
pep_cores_mapping_five['end'] = pep_cores_mapping_five['end'].str.split(';')
pep_cores_mapping_five['accessions'] = pep_cores_mapping_five['accessions'].str.split(';')
pep_cores_mapping_five['proteome_occurrence'] = pep_cores_mapping_five['proteome_occurrence'].str.split(';')
pep_cores_mapping_five['entire_epitope_sequence'] = pep_cores_mapping_five['entire_epitope_sequence'].str.split(';')
pep_cores_mapping_five['core_epitope_sequence'] = pep_cores_mapping_five['core_epitope_sequence'].str.split(';')
pep_cores_mapping_five = pep_cores_mapping_five.explode(['start', 'end', 'accessions', 'proteome_occurrence', 'entire_epitope_sequence', 'core_epitope_sequence']).sort_values(['accessions','start','end','sequence'])
out_cols = list(pep_cores_mapping_five.columns.values)


# read in expected result and sort it
result_large = pd.read_csv(large_result, sep='\t')
result_large = result_large.sort_values(by='sequence').reset_index(drop=True).astype(str)
result_large['sample'] = 'U52'
result_large = result_large[out_cols]
result_large['start'] = result_large['start'].str.split(';')
result_large['end'] = result_large['end'].str.split(';')
result_large['accessions'] = result_large['accessions'].str.split(';')
result_large['proteome_occurrence'] = result_large['proteome_occurrence'].str.split(';')
result_large['entire_epitope_sequence'] = result_large['entire_epitope_sequence'].str.split(';')
result_large['core_epitope_sequence'] = result_large['core_epitope_sequence'].str.split(';')
result_large = result_large.explode(['start', 'end', 'accessions', 'proteome_occurrence', 'entire_epitope_sequence', 'core_epitope_sequence']).sort_values(['accessions','start','end','sequence'])
result_large = result_large.drop_duplicates()

# test if epicore produces expected final result
def test_five():
  assert pep_cores_mapping_five.equals(result_large)