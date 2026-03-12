import os
import pandas as pd
import ast
import click
import logging

from . import __version__
from epicore_utils.modules.compute_cores import compute_consensus_epitopes
from epicore_utils.modules.map_result import map_pep_core, gen_epitope_df
from epicore_utils.modules.visualize_protein import (
    plot_protein_landscape,
    plot_peptide_length_dist,
    plot_core_mapping_peptides_hist,
    qc_plots,
    create_html,
)
from epicore_utils.modules.parse_input import parse_input, proteome_to_dict
from epicore_utils.modules.generate_report import gen_report

logger = logging.getLogger(__name__)


class InputParameter(object):
    """This class contains parameters necessary for the epicore script.

    Attributes:
        min_epi_length (int, optional): An integer of the minimal length of a
            consensus epitope.
        min_overlap (int, optional): An integer of the minimal overlap between
            two epitopes to be grouped to the same consensus epitope.
        max_step_size (int, optional): An integer of the maximal distance
            between the start position of two epitopes to be grouped to the same
            consensus epitope.
        seq_column (str, optional): The string of the header of the column
            containing peptide sequence information in the evidence file.
        protacc_column (str, optional): The string of the header of the column
            containing protein accession information in the evidence file.
        intensity_column (str, optional): The string of the header of the
            column containing intensity information in the evidence file.
        delimiter (str, optional): The delimiter that separates multiple
            entries in one column in the evidence file.
        mod_patter (str, optional): A comma separated string with delimiters
            for peptide modifications
        out_dir (str): A string of the directory, were all output files will be
            saved.
        prot_accession (str, optional): A comma separated string containing the
            protein accession, for which the protein landscape will be visualized.
        start_column (str, optional): The string of the header of the column
            containing the start positions of peptides in proteins.
        end_column (str, optional): The string of the header of the column
            containing the end position of peptides in proteins.
        html (bool, optional): Boolean indicating if a html version of the
            plots should be generated.
        proteome_dict (dict): The specified fasta loaded in a dictionary.
        reference_proteome (str): Path to the specified proteome.
        sample_column (str): The header of the sample column.
        strict (bool, optional): Boolean indicating if the strict mode should be
            run.
        condition_column (str): The header of the column containing condition
            information.
        mapping (bool, optional): Boolean indicating if a mapping of the
            peptides to the computed peptide groups is done.
        included (bool, optional): Boolean indicating if all peptides included
            in the protein region of a peptide group should be added to it.

    """

    def __init__(
        self,
        reference_proteome=None,
        min_epi_length=None,
        min_overlap=None,
        max_step_size=None,
        seq_column=None,
        protacc_column=None,
        intensity_column=None,
        delimiter=None,
        mod_pattern=None,
        out_dir=None,
        prot_accession=None,
        start_column=None,
        end_column=None,
        report=None,
        html=None,
        sample_column=None,
        strict=None,
        condition=None,
        mapping=None,
        included=None,
        qc=None,
        max_group_len=None,
    ):
        self.min_epi_length = min_epi_length
        self.min_overlap = min_overlap
        self.max_step_size = max_step_size
        self.seq_column = seq_column
        self.protacc_column = protacc_column
        self.intensity_column = intensity_column
        self.delimiter = delimiter
        self.mod_pattern = mod_pattern
        self.out_dir = out_dir
        self.prot_accession = prot_accession
        self.start_column = start_column
        self.end_column = end_column
        self.report = report
        self.html = html
        self.proteome_dict = proteome_to_dict(reference_proteome)
        self.reference_proteome = reference_proteome
        self.sample_column = sample_column
        self.strict = strict
        self.condition_column = condition
        self.mapping = mapping
        self.included = included
        self.qc = qc
        self.max_group_len = max_group_len


@click.version_option(__version__, "--version", "-V")
@click.group()
@click.option("--reference_proteome", type=click.Path(exists=True), required=True)
@click.option("--out_dir", type=click.Path(), required=True)
@click.pass_context
def main(ctx, reference_proteome, out_dir):
    ctx.obj = InputParameter(reference_proteome=reference_proteome, out_dir=out_dir)


@click.option("--min_epi_length", type=click.IntRange(5,20), default=11)
@click.option("--min_overlap", type=click.IntRange(1,100), default=9)
@click.option("--max_step_size", type=click.IntRange(0,100), default=5)
@click.option("--seq_column", type=click.STRING, required=True)
@click.option("--sample_column", type=click.STRING, required=True)
@click.option("--protacc_column", type=click.STRING, required=True)
@click.option("--condition_column", type=click.STRING, required=True)
@click.option("--intensity_column", type=click.STRING)
@click.option("--delimiter", type=click.STRING, required=True)
@click.option("--mod_pattern", type=click.STRING)
@click.option("--prot_accession", type=click.STRING)
@click.option("--start_column", type=click.STRING)
@click.option("--end_column", type=click.STRING)
@click.option("--max_group_len", type=click.IntRange(1,100), default=100)
@click.option("--report", is_flag=True)
@click.option("--html", is_flag=True)
@click.option("--strict", is_flag=True)
@click.option("--mapping", is_flag=True)
@click.option("--included", is_flag=True)
@click.option("--QC", is_flag=True)
@click.command()
@click.option("--evidence_file", type=click.Path(exists=True), required=True)
@click.pass_context
def generate_epicore_csv(
    ctx,
    evidence_file,
    min_epi_length,
    min_overlap,
    max_step_size,
    seq_column,
    protacc_column,
    intensity_column,
    delimiter,
    mod_pattern,
    prot_accession,
    start_column,
    end_column,
    report,
    html,
    sample_column,
    strict,
    condition_column,
    mapping,
    included,
    qc,
    max_group_len,
):
    ctx.obj = InputParameter(
        ctx.obj.reference_proteome,
        min_epi_length,
        min_overlap,
        max_step_size,
        seq_column,
        protacc_column,
        intensity_column,
        delimiter,
        mod_pattern,
        ctx.obj.out_dir,
        prot_accession,
        start_column,
        end_column,
        report,
        html,
        sample_column,
        strict,
        condition_column,
        mapping,
        included,
        qc,
        max_group_len,
    )
    if not os.path.exists(ctx.obj.out_dir):
        os.mkdir(ctx.obj.out_dir)
    logging.basicConfig(
        filename=f"{ctx.obj.out_dir}/epicore.log",
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    logger.info(f"Starting epicore run (version: {__version__}) ")
    logger.info(
        f"Parameter: min_epi_length:{min_epi_length}, min_overlap:{min_overlap}, max_step_size:{max_step_size}, max_group_len:{max_group_len}, strict:{strict}, evidence_file:{evidence_file}, fasta_file:{ctx.obj.reference_proteome}, included:{included}, QC:{qc}"
    )

    # ----------------------
    #    Parse input file
    # ----------------------
    protein_df, n_removed_peps, total_intens = parse_input(
        evidence_file,
        ctx.obj.seq_column,
        ctx.obj.protacc_column,
        ctx.obj.intensity_column,
        ctx.obj.start_column,
        ctx.obj.end_column,
        ctx.obj.delimiter,
        ctx.obj.proteome_dict,
        ctx.obj.mod_pattern,
        ctx.obj.sample_column,
        ctx.obj.condition_column,
    )
    protein_df = protein_df.to_pandas()
    os.makedirs(ctx.obj.out_dir, exist_ok=True)

    # ----------------------
    # compute core epitopes, map peptides to cores
    # ----------------------
    protein_df = compute_consensus_epitopes(
        protein_df,
        ctx.obj.min_overlap,
        ctx.obj.max_step_size,
        ctx.obj.min_epi_length,
        ctx.obj.intensity_column,
        ctx.obj.mod_pattern,
        ctx.obj.proteome_dict,
        total_intens,
        ctx.obj.strict,
        ctx.obj.included,
        ctx.obj.max_group_len,
    )
    protein_df.to_csv(f"{ctx.obj.out_dir}/epicore_result.csv")

    if mapping:
        pep_cores_mapping = map_pep_core(
            evidence_file,
            protein_df,
            ctx.obj.seq_column,
            ctx.obj.protacc_column,
            ctx.obj.start_column,
            ctx.obj.end_column,
            ctx.obj.intensity_column,
            ctx.obj.delimiter,
            ctx.obj.mod_pattern,
            ctx.obj.proteome_dict,
        )
        pep_cores_mapping.to_csv(
            f"{ctx.obj.out_dir}/pep_cores_mapping.tsv", sep="\t", index=False
        )

    # ----------------------
    # Reformat data and generate multiple plots
    # ----------------------
    # generate file with one epitope in each row
    epitope_df = gen_epitope_df(protein_df)
    epitope_df.to_csv(f"{ctx.obj.out_dir}/epitopes.csv")

    if qc:
        # plot intern vs extern ratio and consensus sequence coverage
        qc_plots(protein_df, f"{ctx.obj.out_dir}")

    ext = os.path.splitext(evidence_file)[1]
    if ext == ".csv":
        evidence_df = pd.read_csv(evidence_file, delimiter=",", usecols=[ctx.obj.seq_column, ctx.obj.protacc_column])
    elif ext == ".tsv":
        evidence_df = pd.read_csv(evidence_file, delimiter="\t", usecols=[ctx.obj.seq_column, ctx.obj.protacc_column])
    elif ext == ".xlsx":
        evidence_df = pd.read_excel(evidence_file, usecols=[ctx.obj.seq_column, ctx.obj.protacc_column])
    evidence_df[ctx.obj.protacc_column] = evidence_df[ctx.obj.protacc_column].apply(
        lambda accessions: accessions.split(ctx.obj.delimiter)
    )

    # plot number of peptides per peptide group
    fig = plot_core_mapping_peptides_hist(epitope_df)
    fig.savefig(f"{ctx.obj.out_dir}/epitope_intensity_hist.svg")
    if ctx.obj.html:
        create_html(f"{ctx.obj.out_dir}/epitope_intensity_hist.html")

    # plot length distribution of peptides and consensus sequences
    fig, peps, epitopes = plot_peptide_length_dist(
        evidence_df,
        epitope_df,
        ctx.obj.seq_column,
        "consensus_epitopes",
        ctx.obj.seq_column,
        "consensus_epitopes",
        "peptides",
        "consensus sequences",
        mod_pattern,
    )
    fig.savefig(f"{ctx.obj.out_dir}/length_distributions.svg")
    if ctx.obj.html:
        create_html(f"{ctx.obj.out_dir}/length_distributions.html")

        # summarize some results
        if ctx.obj.report:
            gen_report(
                f"http://localhost:8000/{ctx.obj.out_dir}/length_distributions.svg",
                f"http://localhost:8000/{ctx.obj.out_dir}/epitope_intensity_hist.svg",
                epitope_df,
                peps,
                epitopes,
                n_removed_peps,
                ctx,
                evidence_file,
                f"{ctx.obj.out_dir}/epicore_result.csv",
            )


@click.command()
@click.option("--epicore_csv", type=click.Path(exists=True), required=True)
@click.option("--protacc", type=click.STRING, required=True)
@click.pass_context
def plot_landscape(ctx, epicore_csv, protacc):
    if not protacc:
        raise Exception(
            "No protein accession was provided. Please provide a protein accession"
        )
    for accession in protacc.split(","):

        # read in precomputed protein coverage and epitope cores.
        protein_df = pd.read_csv(epicore_csv)

        protein_df["grouped_peptides_start"] = protein_df[
            "grouped_peptides_start"
        ].apply(ast.literal_eval)
        protein_df["core_epitopes_start"] = protein_df["core_epitopes_start"].apply(
            ast.literal_eval
        )
        protein_df["core_epitopes_end"] = protein_df["core_epitopes_end"].apply(
            ast.literal_eval
        )
        protein_df["landscape"] = protein_df["landscape"].apply(ast.literal_eval)

        if accession is not None:
            fig = plot_protein_landscape(protein_df, accession, ctx.obj.proteome_dict)
            fig.savefig(f"{ctx.obj.out_dir}/{accession}.pdf", bbox_inches="tight")
            fig.savefig(f"{ctx.obj.out_dir}/{accession}.svg", bbox_inches="tight")


main.add_command(generate_epicore_csv)
main.add_command(plot_landscape)

if __name__ == "__main__":
    main()
