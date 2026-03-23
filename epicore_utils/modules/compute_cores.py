"""
Computes consensus sequences, by grouping overlapping peptides, building a
landscape for each group and identifying plateaus with a defined minimal length
in each landscape.
"""

import numpy as np
import pandas as pd
import multiprocessing as mp
from multiprocessing import cpu_count


def included_lookback(
    current_group: dict[str, list],
    group_start: list[int],
    group_end: list[int],
    row: pd.Series,
) -> tuple[pd.Series, list[bool]]:
    """Adds peptides to groups in which they are completely included.

    Args:
        current_group: Dictionary containing information about the peptide group.
        group_start: The start positions of the previous peptide groups.
        group_end: The end positions of the previous peptide groups.
        row: The protein and it's peptides that are currently processed.

    Returns:
        A tuple containing the input rows, where peptides are added to the
        peptides groups they are included in, and a boolean array indicating if
        the current peptide is completely included in previous groups.
    """
    included_previous = False
    # iterate over all previous groups which can include the current peptide
    for pos_index, (min_start, max_end) in enumerate(zip(group_start, group_end)):
        if min_start + 23 < int(current_group["starts"][-1]):
            # stop lookback at peptide groups with N-terminal distance > 17 aas
            continue
        else:

            if int(current_group["ends"][-1]) <= max_end:
                # add peptide to previous group
                row["grouped_peptides_start"][-(pos_index + 1)] = row[
                    "grouped_peptides_start"
                ][-(pos_index + 1)] + [int(current_group["starts"][-1])]
                row["grouped_peptides_end"][-(pos_index + 1)] = row[
                    "grouped_peptides_end"
                ][-(pos_index + 1)] + [int(current_group["ends"][-1])]
                row["grouped_peptides_sequence"][-(pos_index + 1)] = row[
                    "grouped_peptides_sequence"
                ][-(pos_index + 1)] + [current_group["sequences"][-1]]
                row["grouped_peptides_sample"][-(pos_index + 1)] = row[
                    "grouped_peptides_sample"
                ][-(pos_index + 1)] + [current_group["samples"][-1]]
                row["grouped_peptides_condition"][-(pos_index + 1)] = row[
                    "grouped_peptides_condition"
                ][-(pos_index + 1)] + [current_group["conditions"][-1]]
                row["peptide_indices"][-(pos_index + 1)] = row["peptide_indices"][
                    -(pos_index + 1)
                ] + [current_group["indices"][-1]]
                included_previous = True

    current_group["included_previous"].append(included_previous)
    return row, current_group


def function_apply(
    df: pd.DataFrame, function: callable, function_args=[]
) -> pd.DataFrame:
    """Applies a function to a dataframe.

    Args:
        df: A dataframe.
        function: A function which should be applied to a dataframe.
        function_args: List of all parameters required by the function.

    Returns:
        The input dataframe to which function was applied to.
    """
    df = df.apply(lambda row: function(row, *function_args), axis=1)
    return df


def parallelized_apply(
    chunk_function: callable, df: pd.DataFrame, function_args=[]
) -> pd.DataFrame:
    """Apply a function to a dataframe using multiprocessing.

    Args:
        chunk_function: The function that should be applied to the chunks.
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
    blocksize = max(len(df) // n_parallel, 1)

    with mp.Pool(min(n_parallel, blocksize)) as pool:
        chunk_dfs = pool.starmap(
            function_apply,
            [
                (
                    df.iloc[
                        chunk
                        * blocksize : (
                            (chunk + 1) * blocksize
                            if chunk < min(n_parallel, blocksize) - 1
                            else len(df)
                        )
                    ],
                    chunk_function,
                    function_args,
                )
                for chunk in range(min(n_parallel, blocksize))
            ],
        )
    df = pd.concat(chunk_dfs)

    if n_proteins != len(df):
        raise Exception("Something went wrong in the multiprocessing. Some rows got \
                        lost. ")

    return df


def group_peptides_protein(
    row: list,
    min_overlap: int,
    max_step_size: int,
    intensity_column: str,
    total_intens: float,
    strict: bool,
    included: bool,
    max_group_len: int,
) -> pd.DataFrame:
    """Group the peptides to consensus epitopes.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        min_overlap: An integer of the minimal overlap between two epitopes
            to be grouped to the same consensus epitope.
        max_step_size: An integer of the maximal distance between the start
            position of two epitopes to be grouped to the same consensus
            epitope.
        intensity_column: The header of the column containing the intensities
            of the peptides.
        total_intens: The total intensity of the evidence file.
        strict: Indicates if strict version should be run.
        included: Indicates if a peptide being included in the protein region of
            a peptide group should be added to it.
        max_group_len: Maximal length of a peptide group.

    Returns:
        The protein_df with the five to seven additional columns:
        grouped_peptides_start, grouped_peptides_end,
        grouped_peptides_sequence, grouped peptides_intensity,
        core_epitopes_intensity, relative_core_intensity and
        sequence_group_mapping. The first five column contain the start
        position, end position, sequence and intensity of the peptides grouped
        to consensus epitopes. The column relative_core_intensity contains the
        relative core intensities of the consensus core epitopes. The column
        sequence group mapping contains for each peptide to which consensus
        epitope it is grouped.
    """

    # Protein data
    start_pos = row["start"]
    end_pos = row["end"]
    sequences = row["sequence"]
    samples = row["sample"]
    conditions = row["condition"]
    if intensity_column:
        intensity = row["intensity"]
    peptide_indices = row["peptide_index"]

    # Stores current peptide group
    current_group = {
        "starts": [],
        "ends": [],
        "sequences": [],
        "samples": [],
        "conditions": [],
        "indices": [],
        "included_previous": [],
        "max_end": 0,
        "prev_end": 0,
        "min_end": 100_000
    }



    if intensity_column:
        grouped_peptides_intensity = []
        core_intensity = 0
    n_jumps = 0
    mapping = []
    group_start = []
    group_end = []

    # Iterate over all peptides mapped to a protein
    for i in range(len(start_pos) - 1):

        # add peptide to current group
        current_group["starts"].append(int(start_pos[i]))
        current_group["ends"].append(int(end_pos[i]))
        current_group["sequences"].append(sequences[i])
        current_group["samples"].append(samples[i])
        current_group["conditions"].append(conditions[i])
        current_group["max_end"] = max(current_group["max_end"], int(end_pos[i]))
        if len(sequences[i]) >= min_overlap: 
            current_group['prev_end'] = int(end_pos[i])
        current_group["indices"].append(peptide_indices[i])
        current_sequence_length = len(sequences[i])
        next_sequence_length = len(sequences[i+1])

        # redefine end position
        if current_sequence_length >= min_overlap: # peptide long enough
            current_group['min_end'] = min(current_group['min_end'], int(end_pos[i]))
        if current_group['min_end'] == 100_000: # enforce new group 
            current_group['min_end'] = 0

        if included:
            # check if peptide is included in previous peptide groups
            row, current_group = included_lookback(
                current_group, group_start, group_end, row
            )

        # check if a new peptide group should be started
        step_size = int(start_pos[i + 1]) - int(start_pos[i])
        mapping.append(n_jumps)
        if intensity_column:
            core_intensity += float(intensity[i])

        overlap = int(end_pos[i]) - int(start_pos[i + 1]) + 1
        if len(sequences[i+1]) < min_overlap:
            if end_pos[i] >= end_pos[i+1]:
                overlap = min_overlap
        if (len(sequences[i]) < min_overlap):
            overlap = current_group['prev_end'] - int(start_pos[i + 1]) + 1

        group_overlap = current_group['min_end'] - int(start_pos[i + 1]) + 1

        # prevent group breakage for sequences shorter then min_overlap
        if (next_sequence_length < min_overlap) and (group_overlap == next_sequence_length):
            group_overlap = min_overlap

        max_len_condition = (
            max_group_len < int(end_pos[i + 1]) - current_group["starts"][0]
        )

        if included:
            condition_group = ((group_overlap < min_overlap) or max_len_condition) and (
                current_group["max_end"] < int(end_pos[i + 1])
            )
        elif strict:
            condition_group = (group_overlap < min_overlap) or max_len_condition
        else:  # loose mode
            condition_group = (
                (step_size >= max_step_size) and (overlap < min_overlap)
            ) or max_len_condition

        # start new peptide group if condition are met
        if step_size != 0:
            if (condition_group and not included) or (
                condition_group
                and included
                and (not all(current_group["included_previous"]))
            ):
                # add group if it is not included in previous group or if included flag is not set
                row["grouped_peptides_start"].append(current_group["starts"])
                row["grouped_peptides_end"].append(current_group["ends"])
                row["grouped_peptides_sequence"].append(current_group["sequences"])
                row["grouped_peptides_sample"].append(current_group["samples"])
                row["grouped_peptides_condition"].append(current_group["conditions"])
                row["peptide_indices"].append(current_group["indices"])
                if intensity_column:
                    row["grouped_peptides_intensity"].append(grouped_peptides_intensity)
                    row["core_epitopes_intensity"].append(core_intensity)

                group_start = [current_group["starts"][0]] + group_start
                group_end = [max(current_group["ends"])] + group_end
            if condition_group:
                # rest current group
                n_jumps += 1
                current_group = {
                    "starts": [],
                    "ends": [],
                    "sequences": [],
                    "samples": [],
                    "conditions": [],
                    "indices": [],
                    "included_previous": [],
                    "max_end": 0,
                    "prev_end": 0,
                    "min_end": 100_000,
                }
                if intensity_column:
                    core_intensity = 0

    # add last peptide to current group
    current_group["starts"].append(int(start_pos[-1]))
    current_group["ends"].append(int(end_pos[-1]))
    current_group["sequences"].append(sequences[-1])
    current_group["samples"].append(samples[-1])
    current_group["conditions"].append(conditions[-1])
    current_group["max_end"] = max(current_group["max_end"], int(end_pos[-1]))
    current_group["indices"].append(peptide_indices[-1])

    # special case for last peptide match of protein
    if len(current_group["ends"]) == 0:

        if included:
            # check if peptide is included in previous peptide groups
            row, current_group = included_lookback(
                current_group, group_start, group_end, row
            )

        if (not included) or (
            included and (not all(current_group["included_previous"]))
        ):
            # add group if it is not included in previous group or if included flag is not set
            row["grouped_peptides_start"].append([int(start_pos[-1])])
            row["grouped_peptides_end"].append([int(end_pos[-1])])
            row["grouped_peptides_sequence"].append([sequences[-1]])
            row["grouped_peptides_sample"].append([samples[-1]])
            row["grouped_peptides_condition"].append([conditions[-1]])
            row["peptide_indices"].append([peptide_indices[-1]])
            if intensity_column:
                row["grouped_peptides_intensity"].append([intensity[-1]])
                row["core_epitopes_intensity"].append(peptide_indices[-1])
            mapping.append(n_jumps)

    else:

        if included:
            # check if peptide is included in previous peptide groups
            row, current_group = included_lookback(
                current_group, group_start, group_end, row
            )

        if intensity_column:
            grouped_peptides_intensity.append(intensity[-1])
            core_intensity += float(intensity[-1])
            row["grouped_peptides_intensity"].append(grouped_peptides_intensity)
            row["core_epitopes_intensity"].append(core_intensity)

        if (not included) or (
            included and (not all(current_group["included_previous"]))
        ):
            # add group if it is not included in previous group or if included flag is not set
            row["grouped_peptides_start"].append(current_group["starts"])
            row["grouped_peptides_end"].append(current_group["ends"])
            row["grouped_peptides_sequence"].append(current_group["sequences"])
            row["grouped_peptides_sample"].append(current_group["samples"])
            row["grouped_peptides_condition"].append(current_group["conditions"])
            row["peptide_indices"].append(current_group["indices"])
            mapping.append(n_jumps)

    row["sequence_group_mapping"] = mapping

    # update start and end position (add start and end position for multiple occurrences)
    row["start"] = [
        f"{start}" for group in row["grouped_peptides_start"] for start in group
    ]
    row["end"] = [f"{end}" for group in row["grouped_peptides_end"] for end in group]
    row["peptide_index"] = [
        f"{index}" for group in row["peptide_indices"] for index in group
    ]

    return row


def group_peptides(
    protein_df: pd.DataFrame,
    min_overlap: int,
    max_step_size: int,
    intensity_column: str,
    total_intens: float,
    strict: bool,
    included: bool,
    max_group_len: int,
) -> pd.DataFrame:
    """Group the peptides based on their overlap.

    Args:
        protein_df: A dataframe containing one protein per row.
        min_overlap: The minimal overlap required between two peptides for them
            to be grouped together.
        max_step_size: The distance between two peptides up to which the overlap
            between the peptide is not considered for the peptide grouping.
        intensity_column: Header of the column containing intensity information.
        total_intens: The total intensity of all peptides in a file.
        strict: Boolean that indicates if the strict mode should be run.
        included: Boolean that indicates if all peptides completely included in
            the protein region of a peptide group should be added to the group.

    Returns:
        A dataframe where each row contains the peptide groups of a protein.
    """

    # start, end, sequence and intensity of peptides of one group grouped together
    protein_df["grouped_peptides_start"] = [[] for _ in range(len(protein_df))]
    protein_df["grouped_peptides_end"] = [[] for _ in range(len(protein_df))]
    protein_df["grouped_peptides_sequence"] = [[] for _ in range(len(protein_df))]
    protein_df["grouped_peptides_sample"] = [[] for _ in range(len(protein_df))]
    protein_df["grouped_peptides_condition"] = [[] for _ in range(len(protein_df))]
    if intensity_column:
        protein_df["grouped_peptides_intensity"] = [[] for _ in range(len(protein_df))]
        # for each peptide group the total and relative intensity of the entire group
        protein_df["core_epitopes_intensity"] = [[] for _ in range(len(protein_df))]
        protein_df["relative_core_intensity"] = [[] for _ in range(len(protein_df))]
    # for each peptide the index of its group
    protein_df["sequence_group_mapping"] = [[] for _ in range(len(protein_df))]
    protein_df["peptide_indices"] = [[] for _ in range(len(protein_df))]

    protein_df = parallelized_apply(
        group_peptides_protein,
        protein_df,
        function_args=[
            min_overlap,
            max_step_size,
            intensity_column,
            total_intens,
            strict,
            included,
            max_group_len,
        ],
    )

    return protein_df


def pep_landscape(row: pd.Series) -> np.array:
    """Generate landscape of a peptide in a group.

    Args:
        row: A row of a pandas dataframe containing per row one protein,
             peptides mapped to the protein with their start and end position
             and their intensity.

    Returns:
        A list containing the landscape of a protein in a given peptide group.
    """
    row["pep_landscape"] = np.zeros([row["end_max"] - row["start_min"] + 1])
    row["pep_landscape"][
        row["grouped_peptides_start"]
        - row["start_min"] : (row["grouped_peptides_end"] - row["start_min"] + 1)
    ] = 1
    return row


def comp_landscape(
    protein_df: pd.DataFrame, proteome_dict: dict[str, str]
) -> pd.DataFrame:
    """Compute the landscape of all consensus epitope groups.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        proteome_dict: A dictionary containing the reference proteome.

    Returns:
        The protein_df with two additional columns (landscape, whole_epitopes).
        The column landscape contains the landscapes of the consensus epitopes
        group. The height of the landscape of a consensus epitope group at a
        position is the number of peptides that contain that position. The
        column whole_epitopes contains the entire sequence of the consensus
        epitope group.
    """
    # get entire sequence of each group
    protein_df["start_min"] = [min(x) for x in protein_df["grouped_peptides_start"]]
    protein_df["end_max"] = [max(x) for x in protein_df["grouped_peptides_end"]]
    protein_df["whole_epitopes"] = protein_df["accession"].map(proteome_dict)
    protein_df["whole_epitopes"] = [
        whole_epitopes[start : end + 1]
        for start, end, whole_epitopes in zip(
            protein_df["start_min"], protein_df["end_max"], protein_df["whole_epitopes"]
        )
    ]  # protein_df['whole_epitopes'].str.slice(start=protein_df['start_min'],stop=protein_df['end_max']+1)
    protein_df = protein_df.assign(whole_epitopes_all=protein_df.whole_epitopes)
    protein_df["group"] = range(len(protein_df))
    protein_df = protein_df.explode(
        [
            "grouped_peptides_start",
            "grouped_peptides_end",
            "grouped_peptides_sample",
            "grouped_peptides_condition",
            "grouped_peptides_sequence",
        ]
    )

    # compute landscape of each peptide
    protein_df = parallelized_apply(pep_landscape, protein_df)

    cols = protein_df.columns.values
    aggr = {
        "whole_epitopes_all": lambda x: list(x),
        "landscape": "first",
        "grouped_peptides_start": lambda x: list(x),
        "grouped_peptides_end": lambda x: list(x),
        "grouped_peptides_sample": lambda x: list(x),
        "grouped_peptides_condition": lambda x: list(x),
        "grouped_peptides_sequence": lambda x: list(x),
    }
    for col in cols:
        if col not in aggr.keys():
            aggr[col] = "first"

    # compute the landscape of the entire group
    comb_landscapes = (
        protein_df.groupby("group")
        .agg({"pep_landscape": lambda column: np.sum(tuple(column), axis=0)})
        .reset_index()
    )
    comb_landscapes = comb_landscapes.rename(columns={"pep_landscape": "landscape"})

    # add group landscapes to the peptide groups
    protein_df = pd.merge(protein_df, comb_landscapes, on="group")
    protein_df["landscape"] = protein_df["landscape"].apply(
        lambda landscape: landscape.tolist()
    )
    protein_df = protein_df.groupby("group").agg(aggr)
    aggr["landscape"] = lambda x: list(x)
    aggr["whole_epitopes"] = lambda x: list(x)
    aggr["whole_epitopes"] = lambda x: list(x)
    del aggr["accession"]
    aggr["start_min"] = lambda x: list(x)

    # combine all groups of one protein
    protein_df = protein_df.groupby("accession").agg(aggr)
    protein_df["whole_epitopes_all"] = protein_df["whole_epitopes_all"].apply(
        lambda x: np.hstack(x).tolist()
    )
    protein_df = protein_df.drop(["group", "pep_landscape"], axis=1)
    protein_df = protein_df.reset_index()

    return protein_df


def find_minima(row: pd.Series) -> np.array:
    """Compute local minima of a landscape.

    Args:
        row: A pandas Series containing one peptide group.

    Returns:
        The input series with the additional information of the minima location.
    """
    landscape = row["landscape"]

    # find minima
    landscape_f = np.roll(landscape, -1) > landscape
    landscape_b = np.roll(landscape, 1) > landscape

    # ensure minima is not at the edge
    edge = np.full(len(landscape), True)
    edge[0] = False
    edge[-1] = False
    landscape = landscape_f & landscape_b & edge

    row["split_rows"] = landscape
    return row


def new_groups(row: pd.Series) -> list:
    """Generate new peptide groups based on the landscape minima.

    Args:
        row: A pandas series containing one peptide group.

    Returns:
        The pandas series with updated peptide groups based on the minima split
        positions.
    """
    splits = row["split_position"]

    row["new_groups_start"] = []
    row["new_groups_end"] = []
    row["new_groups_sample"] = []
    row["new_groups_condition"] = []
    row["new_group_sequences"] = []

    start_list = []
    end_list = []
    sample_list = []
    condition_list = []
    sequence_list = []

    # if the landscape has no minima
    if len(splits) == 0:
        row["new_groups_start"] = [row["grouped_peptides_start"]]
        row["new_groups_end"] = [row["grouped_peptides_end"]]
        row["new_groups_sample"] = [row["grouped_peptides_sample"]]
        row["new_groups_condition"] = [row["grouped_peptides_condition"]]
        row["new_group_sequences"] = [row["grouped_peptides_sequence"]]
        return row

    for start, end, sample, condition, sequence in zip(
        row["grouped_peptides_start"],
        row["grouped_peptides_end"],
        row["grouped_peptides_sample"],
        row["grouped_peptides_condition"],
        row["grouped_peptides_sequence"],
    ):

        if len(splits) == 0 or splits[0] > start:
            # build peptide group before split
            start_list.append(start)
            end_list.append(end)
            sample_list.append(sample)
            condition_list.append(condition)
            sequence_list.append(sequence)

        elif splits[0] <= start:
            # start new peptide group after minimum
            row["new_groups_start"].append(start_list)
            row["new_groups_end"].append(end_list)
            row["new_groups_sample"].append(sample_list)
            row["new_groups_condition"].append(condition_list)
            row["new_group_sequences"].append(sequence_list)
            sequence_list = [sequence]
            start_list = [start]
            end_list = [end]
            sample_list = [sample]
            condition_list = [condition]
            splits = splits[1:]

    row["new_groups_start"].append(start_list)
    row["new_groups_end"].append(end_list)
    row["new_groups_sample"].append(sample_list)
    row["new_groups_condition"].append(condition_list)
    row["new_group_sequences"].append(sequence_list)

    return row


def group_refinement(protein_df: pd.DataFrame, proteome_dict: dict[str, str]):
    """Split peptide groups at landscape minima.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        proteome_dict: A dictionary containing the reference proteome.

    Returns:
        The pandas dataframe with refined peptide groups.
    """

    # convert protein_df so each peptide group is in one row
    protein_df = protein_df.explode(
        [
            "grouped_peptides_start",
            "grouped_peptides_end",
            "grouped_peptides_sample",
            "grouped_peptides_condition",
            "grouped_peptides_sequence",
            "landscape",
            "start_min",
        ]
    )

    # find landscape minima
    protein_df = parallelized_apply(find_minima, protein_df)
    protein_df["split_position"] = protein_df.apply(
        lambda row: np.where(row["split_rows"])[0] + row["start_min"], axis=1
    )

    # refine peptide groups
    protein_df = parallelized_apply(new_groups, protein_df)
    protein_df = protein_df.explode(
        [
            "new_groups_start",
            "new_groups_end",
            "new_groups_sample",
            "new_groups_condition",
            "new_group_sequences",
        ]
    )
    protein_df = protein_df.drop(
        [
            "grouped_peptides_start",
            "grouped_peptides_end",
            "grouped_peptides_sample",
            "grouped_peptides_condition",
            "grouped_peptides_sequence",
        ],
        axis=1,
    )
    protein_df = protein_df.rename(
        columns={
            "new_groups_start": "grouped_peptides_start",
            "new_groups_end": "grouped_peptides_end",
            "new_groups_sample": "grouped_peptides_sample",
            "new_groups_condition": "grouped_peptides_condition",
            "new_group_sequences": "grouped_peptides_sequence",
        }
    )
    protein_df = protein_df.drop(["landscape"], axis=1)

    # recompute group landscapes
    protein_df = comp_landscape(protein_df, proteome_dict)

    return protein_df


def get_consensus_epitopes(row: list, min_epi_len: int) -> pd.DataFrame:
    """Compute the consensus epitope sequence of each consensus epitope group.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        min_epi_length: An integer of the minimal length of a consensus epitope.

    Returns:
        The protein_df with one additional column, that contains the consensus
        epitope sequence of each consensus epitope group.
    """

    # iterate all peptide groups
    for group, landscape in enumerate(row["landscape"]):

        # build consensus epitopes
        total_counts = np.unique(landscape)
        total_counts[::-1].sort()

        # find total coverage for which consensus epitope is at least min_epi_len long
        for total_count in total_counts:

            Z = landscape < total_count

            # get lengths of peptide sequences with coverage above the current threshold
            seqs_idx = np.where(np.diff(np.hstack(([False], ~Z, [False]))))[0].reshape(
                -1, 2
            )
            ce_start_pos = seqs_idx[np.diff(seqs_idx, axis=1).argmax(), 0]
            current_pep_length = np.diff(seqs_idx, axis=1).max()

            # check if min_epi_length is fulfilled for that sequence
            if current_pep_length >= min_epi_len:

                # get position of epitope in protein sequences
                pep_in_prot_start = ce_start_pos.item()
                pep_in_prot_end = pep_in_prot_start + current_pep_length.item()

                # get consensus epitopes
                whole_epitope_wo_mod = row["whole_epitopes"][group]
                for _ in row["grouped_peptides_sequence"][group]:
                    row["consensus_epitopes_all"].append(
                        whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end]
                    )
                    row["core_epitopes_start_all"].append(
                        (pep_in_prot_start + min(row["grouped_peptides_start"][group]))
                    )
                    row["core_epitopes_end_all"].append(
                        pep_in_prot_end + min(row["grouped_peptides_start"][group]) - 1
                    )
                row["consensus_epitopes"].append(
                    whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end]
                )
                row["core_epitopes_start"].append(
                    pep_in_prot_start + min(row["grouped_peptides_start"][group])
                )
                row["core_epitopes_end"].append(
                    pep_in_prot_end + min(row["grouped_peptides_start"][group]) - 1
                )
                break

            # no consensus sequence fulfilling minimal epitope length
            if total_count == total_counts[-1]:
                pep_in_prot_start = ce_start_pos.item()
                pep_in_prot_end = pep_in_prot_start + current_pep_length.item()
                whole_epitope_wo_mod = row["whole_epitopes"][group]
                for _ in row["grouped_peptides_sequence"][group]:
                    row["consensus_epitopes_all"].append(
                        whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end]
                    )
                    row["core_epitopes_start_all"].append(
                        pep_in_prot_start + min(row["grouped_peptides_start"][group])
                    )
                    row["core_epitopes_end_all"].append(
                        pep_in_prot_end + min(row["grouped_peptides_start"][group]) - 1
                    )
                row["consensus_epitopes"].append(
                    whole_epitope_wo_mod[pep_in_prot_start:pep_in_prot_end]
                )
                row["core_epitopes_start"].append(
                    pep_in_prot_start + min(row["grouped_peptides_start"][group])
                )
                row["core_epitopes_end"].append(
                    pep_in_prot_end + min(row["grouped_peptides_start"][group]) - 1
                )

    row["proteome_occurrence"] = [
        row["accession"]
        + ":"
        + str(row["consensus_epitopes_all"][i])
        + ":"
        + str(row["core_epitopes_start_all"][i])
        + "-"
        + str(row["core_epitopes_end_all"][i])
        for i in range(len(row["core_epitopes_start_all"]))
    ]
    return row


def get_consensus_epitopes_protein(
    protein_df: pd.DataFrame, min_epi_length: int
) -> pd.DataFrame:
    """Determine the consensus epitopes of all peptide groups.

    Args:
        protein_df: A dataframe containing one protein per row:
        min_epi_length: The minimal length od the consensus epitopes to be
            identified.

    Returns:
        The input dataframe containing the identified consensus epitopes in
        addition to the input data.
    """

    # start, end, sequence and intensity of peptides of one group grouped together
    protein_df["consensus_epitopes"] = [[] for _ in range(len(protein_df))]
    protein_df["consensus_epitopes_all"] = [[] for _ in range(len(protein_df))]
    protein_df["core_epitopes_start"] = [[] for _ in range(len(protein_df))]
    protein_df["core_epitopes_end"] = [[] for _ in range(len(protein_df))]
    protein_df["core_epitopes_start_all"] = [[] for _ in range(len(protein_df))]
    protein_df["core_epitopes_end_all"] = [[] for _ in range(len(protein_df))]

    # compute the landscape for each peptide in the group
    protein_df = parallelized_apply(
        get_consensus_epitopes, protein_df, function_args=[min_epi_length]
    )

    return protein_df


def reorder_peptides(row: pd.Series, intensity_column: str) -> pd.Series:
    """Reorder the peptides mapped to a protein by their position.

    Args:
        row: A row of a pandas dataframe containing per row one protein,
            peptides mapped to the protein, start and end position of the
            peptide in the protein and the intensity of the peptide.
        intensity_column: The header of the column containing the intensities
            of the peptides.
    Returns:
        A reordered version of the input row, where the start positions are sorted in ascending order, the indices of the other columns are reordered in the same pattern.
    """
    if intensity_column:
        lists = list(
            zip(
                row["start"],
                row["end"],
                row["sequence"],
                row["intensity"],
                row["peptide_index"],
                row["sample"],
                row["condition"],
            )
        )
        lists = sorted(lists, key=lambda x: x[5], reverse=True)
        lists = sorted(lists, key=lambda x: int(x[1]), reverse=True)
        lists = sorted(lists, key=lambda x: x[2], reverse=True)
        sorted_lists = sorted(lists, key=lambda x: int(x[0]))
        starts, ends, sequences, intensities, indices, samples, conditions = zip(
            *sorted_lists
        )
        row[
            [
                "start",
                "end",
                "sequence",
                "intensity",
                "peptide_index",
                "sample",
                "condition",
            ]
        ] = (
            list(starts),
            list(ends),
            list(sequences),
            list(intensities),
            list(indices),
            list(samples),
            list(conditions),
        )
    else:
        lists = list(
            zip(
                row["start"],
                row["end"],
                row["sequence"],
                row["peptide_index"],
                row["sample"],
                row["condition"],
            )
        )
        lists = sorted(lists, key=lambda x: x[4], reverse=True)
        lists = sorted(lists, key=lambda x: x[2], reverse=True)
        lists = sorted(lists, key=lambda x: int(x[1]), reverse=True)
        sorted_lists = sorted(lists, key=lambda x: int(x[0]))
        starts, ends, sequences, indices, samples, conditions = zip(*sorted_lists)
        row[["start", "end", "sequence", "peptide_index", "sample", "condition"]] = (
            list(starts),
            list(ends),
            list(sequences),
            list(indices),
            list(samples),
            list(conditions),
        )
    return row


def compute_consensus_epitopes(
    protein_df: pd.DataFrame,
    min_overlap: int,
    max_step_size: int,
    min_epi_len: int,
    intensity_column: float,
    mod_pattern: str,
    proteome_dict: dict[str, str],
    total_intens: float,
    strict: bool,
    included: bool,
    max_group_len: int,
) -> pd.DataFrame:
    """Compute the core and whole sequence of all consensus epitope groups.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        min_overlap: An integer of the minimal overlap between two epitopes
            to be grouped to the same consensus epitope.
        max_step_size: An integer of the maximal distance between the start
            position of two epitopes to be grouped to the same consensus
            epitope.
        min_epi_len: An integer of the minimal length of a consensus epitope.
        intensity_column: The header of the column containing the intensities
            of the peptides.
        mod_pattern: A comma separated string with delimiters for peptide
            modifications
        proteome_dict: A dictionary containing the reference proteome.
        total_intens: The total intensity of the evidence file.
        strict: A boolean indicating if the strict version should be run.
        included: A boolean indicating if all peptides included in the sequence
            of a peptide group should be included in that peptide group.
        max_group_len: Maximal length of a peptide group.

    Returns:
        The protein_df containing for each protein the core and whole sequence
        of each of its consensus epitope groups.
    """
    # group peptides
    protein_df = parallelized_apply(
        reorder_peptides, protein_df, function_args=[intensity_column]
    )
    protein_df = group_peptides(
        protein_df,
        min_overlap,
        max_step_size,
        intensity_column,
        total_intens,
        strict,
        included,
        max_group_len,
    )
    protein_df = protein_df.explode(
        [
            "grouped_peptides_start",
            "grouped_peptides_end",
            "grouped_peptides_sample",
            "grouped_peptides_sequence",
            "grouped_peptides_condition",
        ]
    )

    # refine peptide groups at landscape minima
    protein_df = comp_landscape(protein_df, proteome_dict)
    if (not strict) and (not included):
        protein_df = protein_df[
            [
                "accession",
                "sequence",
                "start",
                "end",
                "peptide_index",
                "sample",
                "condition",
                "grouped_peptides_start",
                "grouped_peptides_end",
                "grouped_peptides_sequence",
                "grouped_peptides_sample",
                "grouped_peptides_condition",
                "sequence_group_mapping",
                "landscape",
                "whole_epitopes",
                "whole_epitopes_all",
                "start_min",
            ]
        ]
        protein_df = group_refinement(protein_df, proteome_dict)
    protein_df = protein_df[
        [
            "accession",
            "sequence",
            "start",
            "end",
            "peptide_index",
            "sample",
            "condition",
            "grouped_peptides_start",
            "grouped_peptides_end",
            "grouped_peptides_sequence",
            "grouped_peptides_sample",
            "grouped_peptides_condition",
            "sequence_group_mapping",
            "landscape",
            "whole_epitopes",
            "whole_epitopes_all",
        ]
    ]

    # compute consensus sequences
    protein_df = get_consensus_epitopes_protein(protein_df, min_epi_len)
    return protein_df
