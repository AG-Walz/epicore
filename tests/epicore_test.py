# import modules
import pandas as pd
from epicore_utils.modules.parse_input import parse_input, proteome_to_dict
from epicore_utils.modules.map_result import map_pep_core
from epicore_utils.modules.compute_cores import compute_consensus_epitopes
import polars as pl



def run_epicore(proteome_path:str, evidence_path: str, seq_column: str, protacc_column: str, intensity_column: str, start_column: str, end_column: str, delimiter: str, mod_pattern: str, min_overlap: int, max_step_size: int, min_epi_length: int, sample_column, strict, condition_column, included=False, max_group_len=100) -> pd.DataFrame:
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
  protein_df, _, total_intens = parse_input(evidence_path, seq_column , protacc_column, intensity_column, start_column, end_column, delimiter, proteome_dict, mod_pattern, sample_column, condition_column)
  protein_df = protein_df.to_pandas()
  protein_df = compute_consensus_epitopes(protein_df, min_overlap, max_step_size, min_epi_length, intensity_column, mod_pattern, proteome_dict, total_intens, strict, included, max_group_len)
  pep_cores_mapping = map_pep_core(evidence_path,protein_df,seq_column,protacc_column,start_column,end_column,intensity_column,delimiter,mod_pattern, proteome_dict)
  pep_cores_mapping = pep_cores_mapping.sort_values(by='sequence').reset_index(drop=True).astype(str)
  return pep_cores_mapping

# evidence files
evidence_path = 'tests/evidence_file.csv'
evidence_path_three = 'tests/evidence_file_cohort.csv'
evidence_path_four = 'tests/evidence_file_modifications.csv'
evidence_path_five = 'tests/evidence_file_five.csv'
evidence_path_seven = 'tests/evidence_file_seven.csv'
large_evidence = 'tests/large_evidence.csv'
min_landscape_evidence = 'tests/minimal_landscape_evidence.csv'
evidence_path_included = 'tests/evidence_file_included.csv'
evidence_path_limlen = 'tests/evidence_file_lengthlim.csv'

# result files
path_result_one = 'tests/result_one.csv'
path_result_strict = 'tests/result_strict.csv'
path_result_two = 'tests/result_two.csv'
path_result_three = 'tests/result_three.csv'
path_result_four = 'tests/results_four.csv'
path_result_five = 'tests/result_five.csv'
path_result_seven = 'tests/result_seven.csv'
large_result = 'tests/large_evidence_result.csv'
min_landscape_result = 'tests/result_min_landscape.csv'
path_result_included='tests/result_included.csv'
path_result_limlen = 'tests/result_lengthlim.csv'

# fasta files
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
condition_column = 'condition'
strict=False


#############################
# Test one (one sample)
#############################
# run epicore on test file one 
pep_cores_mapping_one = run_epicore(proteome_path, evidence_path, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column)

# read in expected result and sort it
result_file_one = pd.read_csv(path_result_one)
result_file_one = result_file_one.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_one():
  assert pep_cores_mapping_one.equals(result_file_one)

#############################
# Test strict (one sample)
#############################
strict = True
min_overlap = 8
pep_cores_mapping_strict = run_epicore(proteome_path, evidence_path, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column)
out_cols = pep_cores_mapping_strict.columns.values

# read in expected result and sort it
result_file_strict = pd.read_csv(path_result_strict)
result_file_strict = result_file_strict[out_cols].sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_strict():
  assert pep_cores_mapping_strict.equals(result_file_strict)

#############################
# Test two (multiple samples)
#############################
min_overlap = 7
pep_cores_mapping_two = run_epicore(proteome_path, evidence_path_three, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column)
out_cols = pep_cores_mapping_two.columns.values

# read in expected result and sort it
result_file_two = pd.read_csv(path_result_three)
print(result_file_two)
result_file_two = result_file_two[out_cols].sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_two():
  assert pep_cores_mapping_two.equals(result_file_two)


#############################
# Test three (multiple samples, multiple charges, modifications)
#############################

min_epi_length = 3

# run epicore on test file one 
pep_cores_mapping_three = run_epicore(proteome_path, evidence_path_four, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column)
out_cols = pep_cores_mapping_three.columns.values

# read in expected result and sort it
result_file_three = pd.read_csv(path_result_four)
result_file_three = result_file_three[out_cols].sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_three():
  assert pep_cores_mapping_three.equals(result_file_three)


#############################
# Test four 
#############################

max_step_size = 5
min_epi_length = 3

# run epicore on test file one 
pep_cores_mapping_four = run_epicore(proteome_path, evidence_path_five, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column)
pep_cores_mapping_four = pep_cores_mapping_four.sort_values(by='sequence')
out_cols = pep_cores_mapping_four.columns.values

# read in expected result and sort it
result_file_four = pd.read_csv(path_result_five)
result_file_four = result_file_four[out_cols].sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_four():
  assert pep_cores_mapping_four.equals(result_file_four)

#############################
# Test five (one sample, large test)
#############################
strict=False
max_step_size = 5
min_epi_length = 13
min_overlap = 11

# run epicore on test file one 
pep_cores_mapping_five = run_epicore(large_fasta, large_evidence, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column)
pep_cores_mapping_five = pep_cores_mapping_five.sort_values(by='sequence')
pep_cores_mapping_five = pep_cores_mapping_five.loc[:, ~pep_cores_mapping_five.columns.str.contains('^Unnamed')]
pep_cores_mapping_five['start'] = pep_cores_mapping_five['start'].str.split(';')
pep_cores_mapping_five['end'] = pep_cores_mapping_five['end'].str.split(';')
pep_cores_mapping_five['accessions'] = pep_cores_mapping_five['accessions'].str.split(';')
pep_cores_mapping_five['proteome_occurrence'] = pep_cores_mapping_five['proteome_occurrence'].str.split(';')
pep_cores_mapping_five['entire_epitope_sequence'] = pep_cores_mapping_five['entire_epitope_sequence'].str.split(';')
pep_cores_mapping_five['consensus_epitope_sequence'] = pep_cores_mapping_five['consensus_epitope_sequence'].str.split(';')
pep_cores_mapping_five = pep_cores_mapping_five.explode(['start', 'end', 'accessions', 'proteome_occurrence', 'entire_epitope_sequence', 'consensus_epitope_sequence']).sort_values(['accessions','start','end','sequence'])
out_cols = list(pep_cores_mapping_five.columns.values)
pep_cores_mapping_five.reset_index(drop=True,inplace=True)

# read in expected result and sort it
result_large = pd.read_csv(large_result)
result_large = result_large.sort_values(by='sequence').reset_index(drop=True).astype(str)
result_large = result_large.loc[:, ~result_large.columns.str.contains('^Unnamed')]
result_large['sample'] = 'UPN52_1'
result_large = result_large[out_cols]
result_large['start'] = result_large['start'].str.split(';')
result_large['end'] = result_large['end'].str.split(';')
result_large['accessions'] = result_large['accessions'].str.split(';')
result_large['proteome_occurrence'] = result_large['proteome_occurrence'].str.split(';')
result_large['entire_epitope_sequence'] = result_large['entire_epitope_sequence'].str.split(';')
result_large['consensus_epitope_sequence'] = result_large['consensus_epitope_sequence'].str.split(';')
result_large = result_large.explode(['start', 'end', 'accessions', 'proteome_occurrence', 'entire_epitope_sequence', 'consensus_epitope_sequence']).sort_values(['accessions','start','end','sequence'])
result_large = result_large.drop_duplicates()
result_large.reset_index(drop=True,inplace=True)


# test if epicore produces expected final result
def test_five():
  assert pep_cores_mapping_five.equals(result_large)

#############################
# Test six (one sample, minimum in landscape)
#############################
max_step_size = 3
min_epi_length = 11
min_overlap = 6

# run epicore on test file one 
pep_cores_mapping_six = run_epicore(proteome_path, min_landscape_evidence, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column)

# read in expected result and sort it
result_file_six = pd.read_csv(min_landscape_result)
result_file_six = result_file_six.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_six():
  assert pep_cores_mapping_six.equals(result_file_six)



#############################
# Test without start and end position
#############################
# define params
seq_column = 'sequence'
protacc_column = 'accessions'
intensity_column = None
sample_column = 'sample'
start_column = None
end_column = None
delimiter = ';'
mod_pattern = ''
min_overlap = 7
max_step_size = 4
min_epi_length = 10
condition_column = 'condition'
strict=False

# run epicore on test file one 
pep_cores_mapping_seven = run_epicore(proteome_path, evidence_path_seven, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column)

# read in expected result and sort it
result_file_seven = pd.read_csv(path_result_seven)
result_file_seven = result_file_one.sort_values(by='sequence').reset_index(drop=True).astype(str)[['sequence','accessions','charge','score','intensity','sample','condition','start','end','entire_epitope_sequence','consensus_epitope_sequence','proteome_occurrence']]

# test if epicore produces expected final result
def test_seven():
  assert pep_cores_mapping_seven.equals(result_file_seven)



####################
# Test included flag
####################
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
condition_column = 'condition'
strict=True
included=True

# run epicore on test file one 
pep_cores_mapping_included = run_epicore(proteome_path, evidence_path_included, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column, included)

# read in expected result and sort it
result_file_included = pd.read_csv(path_result_included)
result_file_included = result_file_included.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_included():
  assert pep_cores_mapping_included.equals(result_file_included)



##############################
# test group length limitation
##############################
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
condition_column = 'condition'
strict=True
included=True
max_group_len=22

# run epicore on test file one 
pep_cores_mapping_limlen = run_epicore(proteome_path, evidence_path_limlen, seq_column, protacc_column, intensity_column, start_column, end_column, delimiter, mod_pattern, min_overlap, max_step_size, min_epi_length, sample_column, strict, condition_column, included, max_group_len)

# read in expected result and sort it
result_file_limlen = pd.read_csv(path_result_limlen)
result_file_limlen = result_file_limlen.sort_values(by='sequence').reset_index(drop=True).astype(str)

# test if epicore produces expected final result
def test_limlen():
  assert pep_cores_mapping_limlen.equals(result_file_limlen)