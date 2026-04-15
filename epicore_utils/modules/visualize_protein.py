"""
This script contains functions for visualizing the epicore results. The
functions can visualize the consensus sequence coverage, the intern versus the
extern ratio, the number of unique peptides contributing to a peptide group, the
length distribution of the consensus sequences versus the length distribution of
the input peptides and the protein landscape.
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np
import re
import ast
from epicore_utils.modules.compute_cores import parallelized_apply
import logging

logger = logging.getLogger(__name__)


def plot_peptide_length_dist(
    first_df: pd.DataFrame,
    second_df: pd.DataFrame,
    first_explode: str,
    second_explode: str,
    first_column: str,
    second_column: str,
    first_label: str,
    second_label: str,
    mod_pattern: str,
) -> tuple[plt.figure, int, int]:
    """Visualize the distribution of lengths of sequences.

    Args:
        first_df: A pandas dataframe.
        second_df: A pandas dataframe.
        first_explode: The header of the column, that is used to explode the
            dataframe.
        second_column: The header of the column, that is used to explode the
            dataframe.
        first_column: The header of the column, that holds the sequences in
            first_df, for which the length distribution will be plotted.
        second_column: The header of the column, which holds the sequences in
            second_df, for which the length distribution will be plotted.
        first_label: Label for the values of first_column in the plot.
        second_label: Label for the values of second_column in the plot.
        mod_pattern: A comma separated string with delimiters for peptide
            modifications

    Returns:
        A tuple including a matplotlib figure and two integers. The figure is a
        histogram visualizing the length distribution of peptides and epitopes.
        The two integers are the number of peptides and the number of epitopes.
    """
    first_long = first_df.explode(first_explode).drop_duplicates(first_explode)
    second_long = second_df.explode(second_explode).drop_duplicates(
        ["whole_epitopes", "consensus_epitopes"]
    )

    logger.info(
        f"{len(first_long)} peptides were reduced to {len(second_long)} consensus sequences."
    )

    # compute a histogram of the sequence lengths both dataframes
    fig, ax = plt.subplots(layout="constrained")
    seq_first = first_long[first_column]
    seq_second = second_long[second_column]

    # remove modifications and determine sequence length
    pattern = r"\(.*?\)"
    peptide_seqs = seq_first.apply(lambda seq: re.sub(pattern, "", seq))
    pattern = r"\[.*?\]"
    peptide_seqs = peptide_seqs.apply(lambda seq: re.sub(pattern, "", seq))
    if mod_pattern:
        pattern = (
            re.escape(mod_pattern.split(",")[0])
            + r".*?"
            + re.escape(mod_pattern.split(",")[1])
        )
        peptide_seqs = peptide_seqs.apply(lambda seq: re.sub(pattern, "", seq))
    first_len = peptide_seqs.map(lambda pep: len(pep)).to_list()
    second_len = seq_second.map(lambda pep: len(pep)).to_list()

    # plot length histogram
    ax.hist(
        first_len, bins=np.arange(5, 50, 1), color="grey", label=first_label, alpha=0.6
    )
    ax.hist(
        second_len, bins=np.arange(5, 50, 1), color="red", label=second_label, alpha=0.6
    )
    ax.legend()
    ax.set_xlabel("length")
    ax.set_ylabel("count")

    return fig, len(first_long), len(second_long)


def plot_protein_landscape(
    protein_df: pd.DataFrame, accession: str, proteome_dict: dict[str, str]
) -> plt.figure:
    """Visualize the landscape of a protein.

    Args:
        protein_df: A pandas dataframe containing one protein per row.
        accession: The string of a protein accession.
        proteome_dict: A dictionary containing the reference proteome.

    Returns:
        A matplotlib bar plot that visualizes the peptide and core epitope
        distribution across the sequence of the protein with the provided
        accession.
    """

    # get protein sequence
    prot_row = protein_df[(protein_df["accession"] == accession)]
    if len(prot_row) == 0:
        raise Exception("The accession {} is not in your input data.".format(accession))
    prot_seq = proteome_dict[accession]
    prot_landscape = [0 for _ in prot_seq]

    # create figure
    fig_width = max(15, round(len(prot_landscape) / 50))
    fig_height = 3
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), layout="constrained")

    max_intens = 0
    # iterate over all peptide groups
    for group, landscape in enumerate(prot_row["landscape"].iloc[0]):

        # assign the groups alternating colors
        if group % 3 == 0:
            color = "red"
        elif group % 3 == 1:
            color = "blue"
        else:
            color = "green"

        # plot group
        group_start = min(prot_row["grouped_peptides_start"].iloc[0][group])
        for idx, position in enumerate(landscape):
            ax.bar(
                group_start + int(idx), position, width=1, alpha=0.4, color=color
            )  # plot landscape
        for pos in range(
            prot_row["core_epitopes_start"].iloc[0][group],
            prot_row["core_epitopes_end"].iloc[0][group] + 1,
        ):
            ax.bar(
                pos, 0.5, width=1, color=color
            )  # plot consensus sequence of peptide group
        max_intens = max(max_intens, max(landscape))

    ybins = max_intens / 10
    ax.yaxis.set_major_locator(MaxNLocator(nbins=max(ybins, 5)))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=max(fig_width / 10, 25)))
    ax.set_title("Landscape and consensus epitopes of the protein {}".format(accession))
    ax.set_xlabel("Position in protein {}".format(accession))
    ax.set_ylabel("Number of aligned peptides")

    return fig


def plot_core_mapping_peptides_hist(epitope_df: pd.DataFrame) -> plt.figure:
    """A histogram of the number of peptides mapped to each consensus sequence.

    Args:
        epitope_df: A dataframe containing one peptide group per row.

    Returns:
        A histogram visualizing the number of peptides per peptide group.
    """
    fig, ax = plt.subplots(layout="constrained")
    n_peps = epitope_df["grouped_peptides_sequence"].apply(
        lambda sequences: len(set(ast.literal_eval(sequences)))
    )
    ax.hist(n_peps, bins=np.arange(1, max(n_peps) + 2, 1))
    ax.set_yscale("log")
    ax.set_xlabel("Number of peptides contributing to a consensus sequence")
    ax.set_ylabel("Count")
    return fig


def included(
    start: int, end: int, overlap_df: pd.DataFrame, included_df: pd.DataFrame
) -> bool:
    """Checks of a peptide is included in all peptides of the peptide group.

    Args:
        start: The start position of the peptide.
        end: The end position of the peptide.
        overlap_df: Dataframe containing information about peptide overlap.
        included_df: Dataframe containing information about peptide inclusions.

    Returns:
        Boolean that indicates if a peptide overlaps that overlaps with all peptides is completely included in any of the peptides.
    """

    # get peptides that overlap with all peptides in the group
    included_df = included_df[overlap_df[overlap_df.all(axis=1)].index]

    # get row of interest and check if any value is true
    return included_df.loc[f"{start}-{end}"].any()


def get_largest_overlap(
    starts: list[int], ends: list[int], pep_start: int, pep_end: int, intern: bool
) -> int:
    """Computes the maximal overlap of a peptide with the peptides of a group.

    Args:
        starts: The start positions of the peptides in a group.
        ends: The end positions of the peptides in a group.
        pep_start: The start position of a peptide.
        pep_end: The end position of a peptide.
        intern: Boolean indicating if overlap is computed within a group.

    Returns:
        The maximal overlap between the peptide and any peptide of the group.
    """
    max_overlap = 0

    # iterate over all peptides in the peptide group.
    for start, end in zip(starts, ends):
        if (start == pep_start) & (end == pep_end):  # peptide occurs in the group
            if intern:
                continue
            else:
                return 0

        if start <= pep_start:  # peptide occurs before the current
            max_overlap = max(
                min(pep_end - pep_start + 1, end - pep_start + 1), max_overlap
            )

        else:  # peptide occurs after the current
            max_overlap = max(min(pep_end - start + 1, end - start + 1), max_overlap)

    return min(max_overlap, pep_end - pep_start + 1)


def get_minimal_overlap(
    starts: list[int], ends: list[int], pep_start: int, pep_end: int
) -> int:
    """Computes the minimal overlap of a peptide with the peptides of a group.

    Args:
        starts: The start positions of the peptide in one group.
        ends: The end positions of the peptide in one group.
        pep_start: The start position of a peptide.
        pep_end: The end position of a peptide.

    Returns:
        The minimal overlap between the peptide and any peptide of the group.
    """

    if len(starts) == 1:
        return 0

    min_overlap = 100
    # iterate over all peptides in the peptide group.
    for start, end in zip(starts, ends):

        if (start == pep_start) & (end == pep_end):  # peptide occurs in the group
            continue

        if (start <= pep_start) & (end <= pep_end):  # peptide occurs before the current
            min_overlap = max(0, min(end - pep_start + 1, min_overlap))

        elif (start > pep_start) & (end > pep_end):  # peptide occurs after the current
            min_overlap = max(0, min(pep_end - start + 1, min_overlap))

        else:
            min_overlap = min(end - start + 1, min_overlap)

    return min(min_overlap, pep_end - pep_start + 1)


def consensus_coverage(
    start: int, end: int, consensus_start: int, consensus_end: int
) -> int:
    """Computes the consensus sequence coverage by a peptide.

    Args:
        start: The start position of a peptide.
        end: The end position of a peptide.
        consensus_start: The start of the consensus sequence.
        consensus_end: The end of the consensus sequence.

    Returns:
        The consensus sequence coverage.
    """
    consensus_len = consensus_end - consensus_start + 1
    if start <= consensus_start:
        return (
            max(0, min(end - consensus_start + 1, consensus_end - consensus_start + 1))
            / consensus_len
        )
    else:
        return max(0, min(consensus_end - start + 1, end - start + 1)) / consensus_len


def build_overlap_df(
    starts: list[int], ends: list[int]
) -> tuple([pd.DataFrame, pd.DataFrame]):
    """Creates two dataframes tracking peptides overlap and inclusion.

    Args:
        starts: The start positions of the peptides in a peptide group.
        ends: The end positions of the peptides in a peptide group.

    Returns:
        Returns a tuple containing two dataframes. The first dataframe contains
        information about if two peptides overlap with each other. The second
        dataframe contains information about if one peptide is completely
        included in the other.
    """
    # create df tracking if peptide is included or overlaps
    column_names = [f"{start}-{end}" for start, end in zip(starts, ends)]
    overlap_df = pd.DataFrame(columns=column_names)
    included_df = pd.DataFrame(columns=column_names)

    # iterate all peptide pairs of the peptide group
    for i, (start_1, end_1) in enumerate(zip(starts, ends)):
        for start_2, end_2 in zip(starts[i + 1 :], ends[i + 1 :]):

            # check if peptides overlap
            overlap = min(end_2, end_1) - start_2 + 1
            if overlap > 0:
                overlap_bool = True
            else:
                overlap_bool = False
            overlap_df.loc[f"{start_2}-{end_2}", f"{start_1}-{end_1}"] = overlap_bool
            overlap_df.loc[f"{start_1}-{end_1}", f"{start_2}-{end_2}"] = overlap_bool

            # check if peptide is included in other peptide
            included_df.loc[f"{start_2}-{end_2}", f"{start_1}-{end_1}"] = end_2 <= end_1

    included_df.loc[f"{starts[0]}-{ends[0]}", f"{starts[0]}-{ends[0]}"] = None
    return overlap_df, included_df


def all_overlap(row: pd.Series) -> pd.DataFrame:
    """Compute sequence overlap between peptide and peptide of previous and next group.

    Args:
        row: Pandas series containing the peptide groups of one protein.

    Returns:
        A pandas dataframe containing information about the overlap of peptides
        within and outside of a peptide group.
    """

    # define dataframe containing the overlap to the previous, current and next peptide group
    pep = pd.DataFrame(
        columns=[
            "previous",
            "intern_max",
            "intern_min",
            "next",
            "len_pep",
            "len_core",
            "sequence",
            "consensus_sequence_coverage",
            "included",
        ]
    )
    k = 0

    # iterate all peptide groups
    for i, (
        group_starts,
        group_ends,
        group_sequences,
        consensus_start,
        consensus_end,
    ) in enumerate(
        zip(
            row["grouped_peptides_start"],
            row["grouped_peptides_end"],
            row["grouped_peptides_sequence"],
            row["core_epitopes_start"],
            row["core_epitopes_end"],
        )
    ):
        pos_df = pd.DataFrame(
            {"starts": group_starts, "ends": group_ends, "sequences": group_sequences}
        )
        pos_df = pos_df.drop_duplicates()
        group_starts = pos_df["starts"].to_list()
        group_ends = pos_df["ends"].to_list()
        group_sequences = pos_df["sequences"].to_list()
        overlap_df, included_df = build_overlap_df(group_starts, group_ends)

        # iterate all peptide of the current group
        for _, (start, end, sequence) in enumerate(
            zip(group_starts, group_ends, group_sequences)
        ):

            # calculate maximal overlap to previous group
            if i > 0:
                pep.loc[k, "previous"] = get_largest_overlap(
                    row["grouped_peptides_start"][i - 1],
                    row["grouped_peptides_end"][i - 1],
                    start,
                    end,
                    False,
                )
            else:
                pep.loc[k, "previous"] = 0

            # calculate maximal overlap to next group
            if i < len(row["grouped_peptides_start"]) - 1:
                pep.loc[k, "next"] = get_largest_overlap(
                    row["grouped_peptides_start"][i + 1],
                    row["grouped_peptides_end"][i + 1],
                    start,
                    end,
                    False,
                )
            else:
                pep.loc[k, "next"] = 0

            # calculate the maximal and minimal overlap within a group
            pep.loc[k, "intern_max"] = get_largest_overlap(
                group_starts, group_ends, start, end, True
            )
            pep.loc[k, "len_pep"] = end - start + 1
            pep.loc[k, "len_core"] = consensus_end - consensus_start + 1
            pep.loc[k, "sequence"] = sequence
            pep.loc[k, "consensus_sequence_coverage"] = consensus_coverage(
                start, end, consensus_start, consensus_end
            )
            pep.loc[k, "included"] = included(start, end, overlap_df, included_df)
            k += 1

    return pep


def qc_plots(protein_df: pd.DataFrame, out_dir: str):
    """Plot the intern vs extern ratio and consensus sequence coverage plot.

    Args:
        protein_df: A dataframe containing one protein and all its peptide
            groups in a row.
        out_dir: The path to the result directory.
    """

    # calculate the peptide overlap with previous, current and next peptide group
    plot_df = parallelized_apply(all_overlap, protein_df)
    plot_df = pd.concat([*plot_df])
    plot_df["intern_ratio"] = plot_df["intern_max"] / plot_df["len_pep"]
    plot_df["extern_ratio"] = plot_df.apply(
        lambda row: min(1, max(row["previous"], row["next"]) / row["len_pep"]), axis=1
    )
    plot_df.to_csv(f"{out_dir}/ov.csv")

    # plot the intern vs the extern ratio
    figure, axis = plt.subplots(1, 1, dpi=150)
    matrix, xedges, yedges = np.histogram2d(
        plot_df["intern_ratio"], plot_df["extern_ratio"], bins=30
    )
    im1 = axis.imshow(
        np.log(np.flip(matrix.T, 0), where=(np.flip(matrix.T, 0)) != 0,out=np.flip(matrix.T, 0)),
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap="Blues",
    )
    axis.set_xlim([0, 1])
    axis.set_ylim([0, 1])
    cbar = plt.colorbar(im1, shrink=0.5)
    cbar.set_label("Peptide count (log)")
    axis.set_xlabel("Intern ratio")
    axis.set_ylabel("Extern ratio")
    figure.savefig(f"{out_dir}/intern_extern.svg")

    # plot consensus sequence coverage
    figure, axis = plt.subplots(1, 1, dpi=150)
    axis.hist(
        plot_df["consensus_sequence_coverage"],
        100,
        color="red",
        label="all consensus sequences",
    )
    axis.set_yscale("log")
    axis.set_ylabel("Number of peptides")
    axis.set_xlabel(f"Consensus sequence coverage")
    axis.grid()
    plt.tight_layout()
    figure.savefig(f"{out_dir}/consensus_sequence_coverage.svg")


def create_html(out_name: str):
    """Create html version of plot.

    Args:
        out_name: Location where the figure should be stored.
    """
    with open(f"{out_name[:-4]}svg", "r") as svg_file:
        svg_content = svg_file.read()
        svg_content = re.sub(r"<\?xml[^>]+\?>", "", svg_content)
        svg_content = re.sub(r"<!DOCTYPE[^>]+>", "", svg_content)
        html = f"<!DOCTYPE html> <html> <body>{svg_content}</body></html>"
        with open(f"{out_name}", "w") as f:
            f.write(html)
