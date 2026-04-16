"""
Reads in the evidence file and reference profile, computes the peptides
positions in the proteome and links a protein accession with all peptides
associated with the protein.
"""

import pandas as pd
import re
from Bio import SeqIO
import os
import itertools
import polars as pl
from multiprocessing import get_context, cpu_count
from typing import Union

import logging

logger = logging.getLogger(__name__)


def read_entire_id_output(
    id_output: str, polars: bool
) -> Union[pd.DataFrame, pl.DataFrame]:
    """Read in entire evidence file.

    Args:
        id_output: Path to the evidence file.
        polars: Indicates if csv should be loaded in polars or pandas.

    Returns:
        The evidence file in a polars or pandas dataframe.

    Raises:
        Exception: If a file with an unsupported file type is specified as the input
    """
    if polars:
        # determine the file type
        ext = os.path.splitext(id_output)[1]
        if ext == ".csv":
            peptides_df = pl.read_csv(
                id_output, separator=",", infer_schema_length=0
            ).with_row_index("peptide_index")
        elif ext == ".tsv":
            peptides_df = pl.read_csv(
                id_output, separator="\t", infer_schema_length=0
            ).with_row_index("peptide_index")
        elif ext == ".xlsx":
            peptides_df = pl.read_excel(id_output).with_row_index("peptide_index")
        else:
            raise Exception("The file type of your evidence file is not supported. \
                            Please use an evidence file that has one of the \
                            following file types: csv, tsv, xlsx")

    else:
        ext = os.path.splitext(id_output)[1]
        if ext == ".csv":
            peptides_df = pd.read_csv(id_output, delimiter=",")
        elif ext == ".tsv":
            peptides_df = pd.read_csv(id_output, delimiter="\t")
        elif ext == ".xlsx":
            peptides_df = pd.read_excel(id_output)
        else:
            raise Exception("The file type of your evidence file is not supported. \
                            Please use an evidence file that has one of the following \
                            file types: csv, tsv, xlsx")

    return peptides_df


def read_id_output(
    id_output: str,
    seq_column: str,
    protacc_column: str,
    intensity_column: str,
    start_column: str,
    end_column: str,
    delimiter: str,
    sample_column: str,
    condition_column: str,
    polars: bool,
) -> Union[pd.DataFrame, pl.DataFrame]:
    """Read in the evidence file.

    Args:
        id_output: The string of the path to the evidence file.
        seq_column: The string of the header of the column containing
            peptide sequence information in the evidence file.
        protacc_column: The string of the header of the column containing
            protein accession information in the evidence file.
        intensity_column: The string of the header of the column containing
            intensity information in the evidence file.
        start_column: The string of the header of the column containing the
            start positions of peptides in proteins.
        end_column: The string of the header of the column containing the end
            position of peptides in proteins.
        delimiter: The delimiter that separates multiple entries in one column
            in the evidence file.
        sample_column: The header of the column containing the sample
            information.
        condition_column: The header of the column containing the condition
            information.

    Returns:
        A polars dataframe containing the columns sequence, protein accession,
        peptides intensity, sample, condition and optional the columns start
        and end of the input evidence file.

    Raises:
        Exception: If the file type of the provided evidence file is not
            supported.
        Exception: If the provided column does not exist in the provided
            evidence file.
        Exception: If a mandatory column header is not provided.
    """
    peptides_df = read_entire_id_output(id_output, True)
    # check that the mandatory headers are provided and all provided column
    # headers are part of the evidence file
    if seq_column not in peptides_df.columns:
        if seq_column:
            raise Exception(
                "The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.".format(
                    seq_column
                )
            )
        else:
            raise Exception("The header for the column containing the peptide \
                            sequences in the evidence file is mandatory. \
                            Please provide the correct header name.")
    if protacc_column not in peptides_df.columns:
        if protacc_column:
            raise Exception(
                "The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.".format(
                    protacc_column
                )
            )
        else:
            raise Exception("The header for the column containing the protein \
                            accessions in the evidence file is mandatory. \
                            Please provide the correct header name.")
    if (intensity_column not in peptides_df.columns) and intensity_column:
        raise Exception(
            "The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.".format(
                intensity_column
            )
        )
    if (start_column not in peptides_df.columns) and start_column:
        raise Exception(
            "The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.".format(
                start_column
            )
        )
    if (end_column not in peptides_df.columns) and end_column:
        raise Exception(
            "The header {} does not exist in the provided \
                            evidence file. Please provide the correct header.".format(
                end_column
            )
        )

    # read in evidence file (start, stop position and intensity only when
    # provided)
    if start_column and end_column:
        if intensity_column:
            peptides_df = peptides_df.select(
                protacc_column,
                seq_column,
                intensity_column,
                "peptide_index",
                start_column,
                end_column,
                sample_column,
                condition_column,
            )
        else:
            peptides_df = peptides_df.select(
                protacc_column,
                seq_column,
                "peptide_index",
                start_column,
                end_column,
                sample_column,
                condition_column,
            )

        # split if peptide occurs multiple times in proteome
        peptides_df = peptides_df.with_columns(pl.col(start_column).cast(pl.String))
        peptides_df = peptides_df.with_columns(pl.col(end_column).cast(pl.String))
        peptides_df = peptides_df.with_columns(
            (pl.col(start_column).str.split(delimiter).alias(start_column))
        )
        peptides_df = peptides_df.with_columns(
            (pl.col(end_column).str.split(delimiter).alias(end_column))
        )
        peptides_df = peptides_df.with_columns(
            (pl.col(protacc_column).str.split(delimiter).alias(protacc_column))
        )

    else:
        if intensity_column:
            peptides_df = peptides_df.select(
                protacc_column,
                seq_column,
                intensity_column,
                "peptide_index",
                sample_column,
                condition_column,
            )
        else:
            peptides_df = peptides_df.select(
                protacc_column,
                seq_column,
                "peptide_index",
                sample_column,
                condition_column,
            )

        # split accessions if peptide occurs multiple times in proteome
        peptides_df = peptides_df.with_columns(
            (pl.col(protacc_column).str.split(delimiter)).alias(protacc_column)
        )

    return peptides_df


def proteome_to_dict(proteome: str) -> dict[str, str]:
    """Read reference proteome into dictionary.

    Args:
        proteome: The string of the path to the reference proteome.

    Returns:
        The reference proteome as a dictionary.
    """
    proteome_dict = {}
    proteome = SeqIO.parse(open(proteome), "fasta")
    for protein in proteome:
        proteome_dict[protein.id] = str(protein.seq)
    return proteome_dict


def get_pos(accessions: list[str], peptide: str, fasta_dict: dict) -> list[list]:
    """Compute positions of peptide in sequence.

    Args:
        accessions: List of protein accessions to which peptide is mapped.
        peptide: Peptide sequence of interest.
        fasta_dict: Dictionary containing reference proteome.

    Returns:
        Lists containing the accessions and start and end position of a peptide in the proteome.
    """

    starts = ""
    ends = ""
    accessions_str = ""
    for accession in set(accessions):
        if accession != "unmapped":
            seq = fasta_dict[accession]
            groups_pos = re.finditer(f"(?=({peptide}))", seq)
            for group_pos in groups_pos:
                pep_start, pep_end = group_pos.span()
                accessions_str += f";{accession}"
                starts += f";{pep_start}"
                ends += f";{pep_end+len(peptide)-1}"
        else:
            accessions_str = f";unmapped"

    return f"{accessions_str[1:]}~{starts[1:]}~{ends[1:]}"


def add_positions(peptides_df: pl.DataFrame, fasta_dict: dict) -> pl.DataFrame:
    """Adds positions of peptides in the proteins.

    Args:
        peptides_df: Polars dataframe containing the peptide sequence and the proteins they map to.
        fasta_dict: Dictionary containing peptides as keys and accessions as values.

    Returns:
        The input dataframe with peptide positions.
    """
    peptides_df = peptides_df.with_columns(
        pl.struct("accessions", "sequence")
        .map_elements(
            lambda x: get_pos(x["accessions"], x["sequence"], fasta_dict),
            return_dtype=pl.String,
        )
        .str.split("~")
        .alias("positions")
    )
    peptides_df = peptides_df.with_columns(
        pl.col("positions").list.get(0).alias("accessions")
    )
    peptides_df = peptides_df.with_columns(
        pl.col("positions").list.get(1).alias("start")
    )
    peptides_df = peptides_df.with_columns(pl.col("positions").list.get(2).alias("end"))
    peptides_df = peptides_df.with_columns(
        pl.col("accessions").str.split(";").alias("accessions")
    )
    peptides_df = peptides_df.with_columns(
        pl.col("start").str.split(";").alias("start")
    )
    peptides_df = peptides_df.with_columns(pl.col("end").str.split(";").alias("end"))
    peptides_df = peptides_df.drop("positions")
    return peptides_df


def group_repetitive(
    start: list[int],
    end: list[int],
    pep: str,
    acc: str,
    idx: list[int],
    sample: list[str],
    condition: list[str],
) -> str:
    """Group peptide occurrences that belong to the same repetitive region.

    Args:
        start: List of all start positions of the peptides in the protein.
        end: List of all end positions of the peptide in the protein.
        pep: The peptide sequence.
        acc: The protein accession.
        idx: List of all peptide indices.
        sample: List of all samples of the peptides in the protein.
        condition: List of all conditions of the peptides in the protein.

    Returns:
        Returns a tuple (starts,ends), where starts is a list containing the
        updated start positions and ends is a list of updated end positions.
        These contain from the input start and end positions all start and end
        positions that are not part of a repetitive region. For the start and
        end positions that are part of repetitive regions the lowest start
        position and highest end position is kept for each repetitive region.
    """

    current = -1
    updated_start = []
    updated_end = []
    updated_idx = []
    updated_conditions = []
    updated_peps = []
    updated_samples = []

    lists = list(zip(start, end, idx))
    lists = sorted(lists, key=lambda x: int(x[0]))
    start, end, idx = zip(*lists)

    group_ends = []
    # add the first occurrences start positions to the start positions
    for i in idx[0]:
        updated_start.append(str(start[0]))

    # iterate all peptides of the protein
    for pep_pos in range(len(start) - 1):
        for i in idx[pep_pos]:
            group_ends.append(end[pep_pos])
        # two start positions are not part of one repetitive region if the next start position is higher than the current end position
        if int(start[pep_pos + 1]) > int(end[pep_pos]) + 1:  # new group
            for i in idx[pep_pos]:
                updated_idx.append(str(i))
                updated_peps.append(pep)
                updated_samples.append(sample[0])
                updated_conditions.append(condition[0])
                current = pep_pos
            for i in idx[pep_pos + 1]:
                updated_start.append(str(start[pep_pos + 1]))
            # add max end of repetitive group
            for _ in group_ends:
                updated_end.append(str(max(group_ends)))
            group_ends = []

        else:  #  add min_start for repetitive group
            for i in idx[pep_pos]:
                updated_idx.append(str(i))
                updated_peps.append(pep)
                updated_samples.append(sample[0])
                updated_conditions.append(condition[0])
            for i in idx[pep_pos + 1]:
                updated_start.append(str(start[current + 1]))
    for i in idx[-1]:
        group_ends.append(end[-1])
    # add the last occurrences end position to the end positions
    for i in idx[-1]:
        updated_idx.append(str(i))
        updated_peps.append(pep)
        updated_samples.append(sample[0])
        updated_conditions.append(condition[0])
    for _ in group_ends:
        updated_end.append(str(max(group_ends)))

    # reduce each occurrence to one
    updated_df = pd.DataFrame(
        {
            "start": updated_start,
            "end": updated_end,
            "peps": updated_peps,
            "idx": updated_idx,
            "sample": updated_samples,
            "cond": updated_conditions,
        }
    )
    updated_df = updated_df.drop_duplicates()
    updated_start = ";".join(updated_df["start"])
    updated_end = ";".join(updated_df["end"])
    updated_peps = ";".join(updated_df["peps"])
    updated_idx = ";".join(updated_df["idx"])
    updated_samples = ";".join(updated_df["sample"])
    updated_conditions = ";".join(updated_df["cond"])

    return f"{updated_start}|{updated_end}|{updated_idx}|{updated_peps}|{updated_samples}|{updated_conditions}"


def parallelized_apply_polars(
    chunk_function: callable, df: pd.DataFrame, function_args=[]
) -> pd.DataFrame:
    """Apply a function to a dataframe using multiprocessing.

    Args:
        chunk_function: The function that should be applyed to the chunks.
        df: The dataframe which is splitted into chunks.
        function_args: A list containing all arguments of chunk_function.

    Returns:
        The input dataframe to which the chunk_function was applied.

    Raises:
        Exception: If data is lost during multiprocessing.
    """
    n_proteins = len(df)

    # calculate number of maximal processes
    n_parallel = max(1, cpu_count() - 5)
    # size of the chunk_dfs
    block_size = max(len(df) // n_parallel, 1)

    with get_context("spawn").Pool(min(n_parallel, block_size)) as pool:
        chunk_dfs = pool.starmap(
            chunk_function,
            [
                (
                    df[
                        chunk
                        * block_size : (
                            (chunk + 1) * block_size
                            if chunk < min(n_parallel, block_size) - 1
                            else len(df)
                        )
                    ],
                    *function_args,
                )
                for chunk in range(min(n_parallel, block_size))
            ],
        )
    df = pl.concat(chunk_dfs)

    if n_proteins != len(df):
        raise Exception("Something went wrong in the multiprocessing. Some rows \
                        got lost. ")

    return df


def group_repetitive_chunk(
    chunk_df: pl.DataFrame, start: int, end: int
) -> pl.DataFrame:
    """Group repetitive peptides.

    Args:
        chunk_df: A polars DataFrame containing one protein per row.
        start: The row, at which the polars DataFrame gets sliced.
        end: The row, at which the polars DataFrame gets sliced.

    Returns:
        A slice of the input DataFrame with grouped repetitive peptides.
    """
    return chunk_df[start:end].with_columns(
        pl.struct(
            "start",
            "end",
            "sequence",
            "accessions",
            "peptide_index",
            "sample",
            "condition",
        )
        .map_elements(
            lambda x: group_repetitive(
                x["start"],
                x["end"],
                x["sequence"],
                x["accessions"],
                x["peptide_index"],
                x["sample"],
                x["condition"],
            ),
            return_dtype=pl.String,
        )
        .str.split("|")
        .alias("repetitive")
    )


def prot_pep_link(
    peptides_df: pd.DataFrame,
    intensity_column: str,
    delimiter: str,
) -> pl.DataFrame:
    """Converts a dataframe from one peptide per row to one protein per row.

    Args:
        peptides_df: A pandas dataframe containing the columns peptide sequence,
            protein accession, intensity, start position and end position.
        intensity_column: The string of the header of the column containing
            intensity information in the evidence file.
        delimiter: The delimiter that separates multiple entries in one column
            in the evidence file.

    Returns:
        A polars dataframe containing one protein per row and all peptides
        mapped to that protein in the peptides_df, with their start position,
        end position and intensity.

    Raises:
        Exception: If a peptide does not occur in the protein to which it is
            mapped.
    """

    if intensity_column:
        proteins = pd.DataFrame(
            columns=[
                "accession",
                "sequence",
                "intensity",
                "start",
                "end",
                "peptide_index",
            ]
        )
    else:
        proteins_df = peptides_df.explode("accessions", "start", "end")
        proteins_df = proteins_df.with_columns(pl.col("start").cast(pl.Int64))

        proteins_df = proteins_df.group_by(
            "sequence", "accessions", "sample", "condition"
        ).agg(pl.col("start"), pl.col("end"), pl.col("peptide_index"))

        proteins_df = proteins_df.with_columns(pl.col("end").cast(pl.List(pl.Int64)))
        proteins_df = proteins_df.with_columns(pl.col("start").cast(pl.List(pl.Int64)))
        proteins_df = proteins_df.with_columns(
            pl.col("peptide_index").cast(pl.List(pl.List(pl.Int64)))
        )
        proteins_df = proteins_df.with_columns(
            pl.col("sample").cast(pl.List(pl.String))
        )
        proteins_df = proteins_df.with_columns(
            pl.col("condition").cast(pl.List(pl.String))
        )

        n_parallel = max(1, cpu_count() - 5)
        block_size = max(len(proteins_df) // n_parallel, 1)
        with get_context("spawn").Pool(n_parallel) as pool:
            chunk_dfs = pool.starmap(
                group_repetitive_chunk,
                [
                    (
                        proteins_df,
                        chunk * block_size,
                        (
                            (chunk + 1) * block_size
                            if chunk < (n_parallel - 1)
                            else len(proteins_df)
                        ),
                    )
                    for chunk in range(n_parallel)
                ],
            )
        proteins_df = pl.concat(chunk_dfs)

        proteins_df = proteins_df.with_columns(
            pl.col("repetitive").list.get(0).alias("start")
        )
        proteins_df = proteins_df.with_columns(
            pl.col("repetitive").list.get(1).alias("end")
        )
        proteins_df = proteins_df.with_columns(
            pl.col("repetitive").list.get(2).alias("peptide_index")
        )
        proteins_df = proteins_df.with_columns(
            pl.col("repetitive").list.get(3).alias("sequence")
        )
        proteins_df = proteins_df.with_columns(
            pl.col("repetitive").list.get(4).alias("sample")
        )
        proteins_df = proteins_df.with_columns(
            pl.col("repetitive").list.get(5).alias("condition")
        )
        proteins_df = proteins_df.with_columns(
            (pl.col("start").str.split(delimiter)).alias("start")
        )
        proteins_df = proteins_df.with_columns(
            (pl.col("end").str.split(delimiter)).alias("end")
        )
        proteins_df = proteins_df.with_columns(
            (pl.col("peptide_index").str.split(delimiter)).alias("peptide_index")
        )
        proteins_df = proteins_df.with_columns(
            (pl.col("sequence").str.split(delimiter)).alias("sequence")
        )
        proteins_df = proteins_df.with_columns(
            (pl.col("sample").str.split(delimiter)).alias("sample")
        )
        proteins_df = proteins_df.with_columns(
            (pl.col("condition").str.split(delimiter)).alias("condition")
        )

        proteins_df = proteins_df.explode(
            "start", "end", "peptide_index", "sequence", "sample", "condition"
        )
        proteins_df = proteins_df.group_by("accessions").agg(
            pl.col("sequence"),
            pl.col("start"),
            pl.col("end"),
            pl.col("peptide_index"),
            pl.col("sample"),
            pl.col("condition"),
        )

        proteins_df = proteins_df.rename({"accessions": "accession"})

    return proteins_df


def parse_input(
    evidence_file: str,
    seq_column: str,
    protacc_column: str,
    intensity_column: str,
    start_column: str,
    end_column: str,
    delimiter: str,
    proteome_dict: dict[str, str],
    mod_pattern: str,
    sample_column: str,
    condition_column: str,
) -> tuple[pd.DataFrame, int, float]:
    """Parse the evidence file.

    Args:
        evidence_file: The string of the path to the evidence file.
        seq_column: The string of the header of the column containing
            peptide sequence information in the evidence file.
        protacc_column: The string of the header of the column containing
            protein accession information in the evidence file.
        intensity_column: The string of the header of the column containing
            intensity information in the evidence file.
        start_column: The string of the header of the column containing the
            start positions of peptides in proteins.
        end_column: The string of the header of the column containing the end
            position of peptides in proteins.
        delimiter: The delimiter that separates multiple entries in one column
            in the evidence file.
        proteome_dict: A dictionary containing the reference proteome.
        mod_pattern: A comma separated string with delimiters for peptide
            modifications

    Returns:
        A tuple of a pandas dataframe, an integer and a float. The dataframe contains
        one protein per row and all peptides mapped to that protein in the
        peptides_df, with their start position, end position and intensity.
        The integer is the number of peptides that are removed due to their
        accession not appearing in the proteome file. The float is the total intensity
        of the evidence file.
    """
    peptides_df = read_id_output(
        evidence_file,
        seq_column,
        protacc_column,
        intensity_column,
        start_column,
        end_column,
        delimiter,
        sample_column,
        condition_column,
        True,
    )
    peptides_df = peptides_df.rename(
        {
            seq_column: "sequence",
            protacc_column: "accessions",
            sample_column: "sample",
            condition_column: "condition",
        }
    )
    if start_column and end_column:
        peptides_df = peptides_df.rename({start_column: "start", end_column: "end"})

    # get peptides/proteins with protein accessions that do not appear in the proteome
    peptides = pl.Series(
        peptides_df.with_columns(
            (
                pl.col("accessions").list.filter(
                    ~pl.element().is_in(list(proteome_dict.keys()))
                )
            ).alias("removed")
        ).select("removed")
    ).to_list()
    n_removed_proteins = set(itertools.chain.from_iterable(peptides))

    # remove all peptides occurring multiple times in different modification and charge states
    peptides_df = peptides_df.with_columns(
        (pl.col("sequence").str.replace_all(r"\(.*?\)", "")).alias("sequence")
    )
    # remove peptides with protein accessions that do not appear in the proteome
    if start_column and end_column:
        peptides_df = peptides_df.with_columns(
            pl.col("start").list.gather(
                pl.col("accessions").list.eval(
                    pl.arg_where(~pl.element().is_in(n_removed_proteins)).alias("start")
                )
            )
        )
        peptides_df = peptides_df.with_columns(
            pl.col("end").list.gather(
                pl.col("accessions").list.eval(
                    pl.arg_where(~pl.element().is_in(n_removed_proteins)).alias("end")
                )
            )
        )
        peptides_df = peptides_df.group_by(
            ["accessions", "sequence", "sample", "condition"]
        ).agg(
            pl.col("start").first(),
            pl.col("end").first(),
            pl.col("peptide_index"),
        )
    else:
        peptides_df = peptides_df.group_by(
            ["accessions", "sequence", "sample", "condition"]
        ).agg(pl.col("peptide_index"))
    peptides_df = peptides_df.with_columns(
        pl.col("accessions").list.gather(
            pl.col("accessions").list.eval(
                pl.arg_where(~pl.element().is_in(n_removed_proteins)).alias(
                    "accessions"
                )
            )
        )
    )
    # remove peptides that are not annotated with any proteome accession
    peptides_df = peptides_df.remove(pl.col("accessions").list.len() == 0)

    # compute start and end positions of the peptides
    if not start_column and not end_column:
        peptides_df = add_positions(peptides_df, proteome_dict)

    logger.info(
        f"Peptides mapped to the following {len(n_removed_proteins)} proteins were removed since the proteins do not appear in the proteome fasta file: {n_removed_proteins}."
    )
    protein_df = prot_pep_link(peptides_df, intensity_column, delimiter)
    if intensity_column:
        total_intens = peptides_df[intensity_column].sum()
    else:
        total_intens = 0

    return protein_df, len(n_removed_proteins), total_intens
