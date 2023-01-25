import re

from pathlib import Path

from cogent3 import get_app, load_table, open_data_store
from mdeq.utils import (
    load_from_sqldb,
    rich_display,
    summary_not_completed,
    write_to_sqldb,
)
from rich.progress import track
from scitrack import CachingLogger

from mdeq_analysis.align_checked import sequence_filter


_suffixes = re.compile(r"(\.fa|\.gz|.json)")


def get_chrom1_ids(cds_indir, intron_indir, metadata_path) -> set:
    """returns set of id's present in all 3 sources"""
    table = load_table(metadata_path)
    table = table.filtered(
        lambda x: "Homo sapiens:chromosome:1" in x, columns="location"
    )
    chrom1_ids = set(table.columns["refid"])

    cds_dstore = open_data_store(cds_indir, suffix=".fa.gz")
    cds_ids = {_suffixes.sub("", m.unique_id) for m in cds_dstore.completed}

    intron_dstore = open_data_store(intron_indir, suffix=".fa.gz")
    intron_ids = {_suffixes.sub("", m.unique_id) for m in intron_dstore.completed}
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
    all_species = get_app("take_named_seqs", "Human", "Chimp", "Gorilla")
    loader = get_app("load_unaligned", moltype="dna", format="fasta")
    aligner = get_app("progressive_align", "codon")

    if outpath.exists() and not overwrite:
        raise IOError(f"{str(outpath)} exists")

    out_dstore = open_data_store(outpath, mode="w")

    writer = write_to_sqldb(out_dstore)
    dstore = open_data_store(cds_indir, suffix=".fa.gz", limit=limit)
    selected = [
        m for m in dstore.completed if _suffixes.sub("", m.unique_id) in selected_ids
    ]
    app = loader + all_species + aligner + writer
    out_dstore = app.apply_to(
        selected, logger=LOGGER, cleanup=True, show_progress=verbose
    )

    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

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
    cds_dstore = open_data_store(cds_path, limit=limit)

    # for processing the CDS sequences
    loader = load_from_sqldb()
    ape_names = "Human", "Chimp", "Gorilla"
    take_seqs = get_app("take_named_seqs", *ape_names)
    seq_qual = sequence_filter()
    third_pos = get_app("take_codon_positions", 3)
    no_degenerates = get_app("omit_degenerates", moltype="dna", gap_is_degen=True)
    no_shorties = get_app("min_length", 300)

    cds_app = loader + take_seqs + seq_qual + third_pos + no_degenerates + no_shorties
    cds_out = outdir / f"filtered-{cds_path.stem}.sqlitedb"
    name = cds_path.stem.replace("cds", "intron")
    intron_out = outdir / f"filtered-{name}.sqlitedb"

    for outpath in (cds_out, intron_out):
        if outpath.exists() and not overwrite:
            raise IOError(f"{str(outpath)} exists")

    cds_out_dstore = open_data_store(cds_out, mode="w")
    cds_writer = write_to_sqldb(cds_out_dstore)

    # for processing the intron sequences, which are in a directory
    intron_dstore = open_data_store(intron_indir, suffix="fa.gz", limit=limit)
    loader = get_app("load_aligned", format="fasta", moltype="dna")
    take_seqs = get_app("take_named_seqs", *ape_names)
    seq_qual = sequence_filter()
    no_degenerates = get_app("omit_degenerates", moltype="dna", gap_is_degen=True)
    no_shorties = get_app("min_length", 3000)

    intron_app = loader + take_seqs + seq_qual + no_degenerates + no_shorties
    intron_out_dstore = open_data_store(intron_out, mode="w")
    intron_writer = write_to_sqldb(intron_out_dstore)
    # the id's have been matched already through the codon-align step
    for cds_m in track(cds_dstore.completed):
        name = _suffixes.sub("", cds_m.unique_id)
        intrn_m = [m for m in intron_dstore if name in m.unique_id]
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
        cds_writer.main(ncomp)
        intron_writer.main(ncomp)

    LOGGER.shutdown()
    log = log_file_path.read_text()
    intron_writer.data_store.write_log(data=log, unique_id=log_file_path.name)
    cds_writer.data_store.write_log(data=log, unique_id=log_file_path.name)

    if verbose:
        rich_display(cds_writer.data_store.describe)
        rich_display(intron_writer.data_store.describe)

    log_file_path.unlink()
    return True
