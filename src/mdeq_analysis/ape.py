import re

from pathlib import Path

from cogent3 import load_table
from cogent3.app import align, io, sample
from mdeq.utils import rich_display
from scitrack import CachingLogger


_suffixes = re.compile(r"(\.fa|\.gz)")


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


