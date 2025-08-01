from pathlib import Path

from cogent3 import get_app, make_unaligned_seqs, open_data_store
from cogent3.app import typing as c3_types
from cogent3.app.composable import define_app
from cogent3.parse.fasta import MinimalFastaParser
from mdeq.utils import (
    load_from_sqldb,
    rich_display,
    summary_not_completed,
    write_to_sqldb,
)
from scitrack import CachingLogger

# we combine the unaligned intron sequences from 3 separate directories
# we also rename the sequences

# samples the sequences from Ensembl
# need function to match seqs from different species / introns, producing a merged collection
# rename sequences
# align the sequences
# attach / transfer annotations


def _rename_rodent_seqs(orig):
    orig = orig.lower()
    if "spret" in orig:
        return "msp"
    elif "musculus" in orig:
        return "mmu"
    elif "rat" in orig:
        return "rno"
    else:
        raise NotImplementedError(f"unexpected name {orig}")


@define_app
def rodent_rename(seq_coll: c3_types.UnalignedSeqsType) -> c3_types.UnalignedSeqsType:
    """renames Mus musculus to mmu, Mus spretus to msp and rattus norvegicus to rno"""
    return seq_coll.rename_seqs(_rename_rodent_seqs) if seq_coll else seq_coll


def group_orthologs(
    base_name: str,
    mmu_path: Path,
    msp_path: Path,
    rno_path: Path,
    logger: CachingLogger,
):
    """groups orthologous sequences from the 3 rodent species into a single SequenceCollection

    Parameters
    ----------
    base_name : str
        ortholog filename common across all directories
    mmu_path : Path
        Mus musculus directory
    msp_path : Path
        Mus spretus directory
    rno_path : Path
        Rattus norvegicus directory
    logger:
        for logging input files

    Returns
    -------
    SequenceCollection
    """
    species_dirs = mmu_path, msp_path, rno_path
    data = []
    for dname in species_dirs:
        if logger:
            logger.input_file(dname / base_name)
        data.extend(MinimalFastaParser(dname / base_name))

    return make_unaligned_seqs(
        data,
        moltype="dna",
        info={"source": base_name},
    )


def make_aligned(mmu_path, msp_path, rno_path, outpath, parallel, overwrite, verbose):
    """aggregates homologs and aligns them

    Parameters
    ----------
    mmu_path : Path
        directory containing fasta files with single sequences for Mouse
    msp_path : Path
        directory containing fasta files with single sequences for Algerian Mouse
    rno_path : Path
        directory containing fasta files with single sequences for the Rat
    outpath : Path
        alignment output sqlitedb file path
    parallel : bool
        run in parallel
    overwrite : bool
        overwites existing results
    verbose : bool
        control information display

    Notes
    -----
    Assumes files of homologs from the different species have the same name.
    Sequences are renamed to mmu, msp andf rno for Mouse, Algerian Mouse and
    Rat
    """
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()
    LOGGER.log_file_path = outpath.parent / "fxy-align.log"

    basenames = [p.name for p in mmu_path.glob("*.fasta")]
    args = (mmu_path, msp_path, rno_path, LOGGER)
    seq_collections = [group_orthologs(bn, *args) for bn in basenames]

    renamer = rodent_rename()
    aligner = get_app(
        "progressive_align",
        "nucleotide",
        guide_tree="((msp:0.01,mmu:0.01):0.02,rno:0.03)",
    )
    out_dstore = open_data_store(outpath, mode="w")
    writer = write_to_sqldb(out_dstore)
    app = renamer + aligner + writer
    out_dstore = app.apply_to(
        seq_collections,
        parallel=parallel,
        show_progress=verbose > 0,
        logger=LOGGER,
        cleanup=True,
    )

    out_dstore.unlock()
    success = len(out_dstore.completed) > 1
    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    return success


def filter_positions(inpath, outpath, overwrite, verbose):
    """filters alignment positions, removing non-canonical bases"""
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    LOGGER.log_file_path = outpath.parent / "mdeqasis-fxy-filter_alignments.log"
    if inpath.suffix == ".sqlitedb":
        dstore = open_data_store(inpath)
        loader = load_from_sqldb()
    else:
        dstore = open_data_store(inpath, suffix="fa")
        loader = get_app("load_aligned", moltype="dna", format="fasta")

    just_nucs = get_app(
        "omit_degenerates", moltype="dna", motif_length=1, gap_is_degen=True
    )

    if outpath.exists() and not overwrite:
        return outpath

    outpath.parent.mkdir(parents=True, exist_ok=True)
    out_dstore = open_data_store(outpath, mode="w")
    writer = write_to_sqldb(out_dstore)

    app = loader + just_nucs + writer
    out_dstore = app.apply_to(
        dstore.completed, cleanup=True, show_progress=True, logger=LOGGER
    )
    out_dstore.unlock()
    success = len(out_dstore.completed) > 0
    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    out_dstore.close()
    return success
