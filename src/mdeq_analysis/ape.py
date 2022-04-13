import re

from pathlib import Path

from cogent3 import load_table
from cogent3.app import align, io, sample
from mdeq.sqlite_data_store import sql_loader, sql_writer
from mdeq.utils import rich_display
from rich.progress import track
from scitrack import CachingLogger

from mdeq_analysis import align_check


_suffixes = re.compile(r"(\.fa|\.gz|.json)")


def get_chrom1_ids(cds_indir, intron_indir, metadata_path) -> set:
    """returns set of id's present in all 3 sources"""
    table = load_table(metadata_path)
    table = table.filtered(
        lambda x: "Homo sapiens:chromosome:1" in x, columns="location"
    )
    chrom1_ids = set(table.columns["refid"])

    cds_dstore = io.get_data_store(cds_indir, suffix=".fa.gz")
    cds_ids = {_suffixes.sub("", m.name) for m in cds_dstore}

    intron_dstore = io.get_data_store(intron_indir, suffix=".fa.gz")
    intron_ids = {_suffixes.sub("", m.name) for m in intron_dstore}
    return chrom1_ids & cds_ids & intron_ids


def codon_align(
    cds_indir, intron_indir, metadata_path, outpath, limit, overwrite, verbose, testrun
):
    """codon aligns only files with matching data in the intron data set and on chrom 1"""
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    outpath = Path(outpath)

    LOGGER.log_file_path = outpath.parent / "mdeqasis-dros_sample_alignments.log"
    LOGGER.input_file(metadata_path)

    selected_ids = get_chrom1_ids(cds_indir, intron_indir, metadata_path)
    all_species = sample.take_named_seqs("Human", "Chimp", "Gorilla")
    loader = io.load_unaligned(moltype="dna", format="fasta")
    aligner = align.progressive_align("codon")
    writer = io.write_db(
        outpath, create=True, if_exists=io.OVERWRITE if overwrite else io.RAISE
    )
    dstore = io.get_data_store(cds_indir, suffix=".fa.gz", limit=limit)
    dstore = dstore.filtered(
        callback=lambda x: _suffixes.sub("", x.name) in selected_ids
    )
    app = loader + all_species + aligner + writer
    _ = app.apply_to(dstore, logger=LOGGER, cleanup=True, show_progress=verbose)
    rich_display(app.data_store.describe)
    if len(app.data_store.incomplete) > 0 and verbose:
        rich_display(app.data_store.summary_incomplete)

    return True


def match_cds_intron(
    cds_path, intron_indir, outdir, limit, overwrite, verbose, testrun
):
    """paired sampling of introns with their cds sequences that satisfy filtering standard

    Parameters
    ----------
    cds_path : Path
        path to directory containing alignments of cds sequences
    intron_indir : Path
        directory containing alignments of intron sequences
    outdir : Path
        writes autonamed output for cds and intron here, prefixes each with
        filtered-
    limit : int
        number to sample from each
    overwrite : bool
        replace existing files
    verbose : bool
        messages to screen
    testrun : bool
        limit sampling
    """
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()
    log_file_path = outdir / f"{cds_path.stem}.log"
    LOGGER.log_file_path = log_file_path
    LOGGER.input_file(cds_path)
    cds_dstore = io.get_data_store(cds_path, limit=limit)

    # for processing the CDS sequences
    loader = sql_loader()
    ape_names = "Human", "Chimp", "Gorilla"
    take_seqs = sample.take_named_seqs(*ape_names)
    align_qual = align_check.alignment_filter()
    seq_qual = align_check.sequence_filter()
    third_pos = sample.take_codon_positions(3)
    no_degenerates = sample.omit_degenerates(moltype="dna", gap_is_degen=True)
    no_shorties = sample.min_length(300)
    cds_app = loader + take_seqs + seq_qual + third_pos + no_degenerates + no_shorties
    cds_out = outdir / f"filtered-{cds_path.stem}.sqlitedb"
    name = cds_path.stem.replace("cds", "intron")
    intron_out = outdir / f"filtered-{name}.sqlitedb"
    cds_writer = sql_writer(
        cds_out,
        create=True,
        if_exists=io.OVERWRITE if overwrite else io.RAISE,
    )

    # for processing the intron sequences, which are in a directory
    intron_dstore = io.get_data_store(intron_indir, suffix="fa.gz", limit=limit)
    loader = io.load_aligned(format="fasta", moltype="dna")
    take_seqs = sample.take_named_seqs(*ape_names)
    align_qual = align_check.alignment_filter()
    seq_qual = align_check.sequence_filter()
    no_degenerates = sample.omit_degenerates(moltype="dna", gap_is_degen=True)
    no_shorties = sample.min_length(3000)
    intron_app = loader + take_seqs + seq_qual + no_degenerates + no_shorties
    intron_writer = sql_writer(
        intron_out,
        create=True,
        if_exists=io.OVERWRITE if overwrite else io.RAISE,
    )
    # the id's have been matched already through the codon-align step
    for cds_m in track(cds_dstore):
        name = _suffixes.sub("", cds_m.name)
        intrn_m = intron_dstore.filtered(pattern=f"*{name}*")
        if not intrn_m:
            continue

        cds = cds_app(cds_m)
        intron = intron_app(intrn_m[0]) if cds else None
        if cds and intron:
            cds_writer(cds)
            intron_writer(intron)
            continue

        # one of these will be a not-completed record, which evaluates to False
        ncomp = intron if cds else cds
        ncomp.message = f"EITHER INTRON OR CDS SAMPLING FAILED\n{ncomp.message}"
        # todo writers should correctly handle writing NotCompleted objects
        # without requiring accessing the data_store attribute
        cds_writer.data_store.write(ncomp)
        intron_writer.data_store.write(ncomp)

    LOGGER.shutdown()
    log = Path(log_file_path).read_text()
    intron_writer.data_store.add_log(log)
    cds_writer.data_store.add_log(log)

    if verbose:
        rich_display(cds_writer.data_store.describe)
        rich_display(intron_writer.data_store.describe)

    return True
