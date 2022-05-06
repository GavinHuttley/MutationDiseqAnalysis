from pathlib import Path

from cogent3 import make_unaligned_seqs
from cogent3.app import io
from cogent3.app.align import progressive_align
from cogent3.app.composable import (
    ALIGNED_TYPE,
    SEQUENCE_TYPE,
    SERIALISABLE_TYPE,
    appify,
)
from cogent3.parse.fasta import MinimalFastaParser
from mdeq.sqlite_data_store import sql_loader, sql_writer
from mdeq.utils import rich_display
from scitrack import CachingLogger


# we combine the unaligned intron sequences from 3 separate directories
# we also rename the sequences

# samples the sequences from Ensembl
# need function to match seqs from different species / introns, producing a merged collection
# rename sequences
# align the sequences
# attach / transfer annotations

_types = (SERIALISABLE_TYPE, SEQUENCE_TYPE, ALIGNED_TYPE)


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


@appify(_types, _types)
def rodent_rename(seq_coll):
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
        data=data,
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
    aligner = progressive_align(
        "nucleotide", guide_tree="((msp:0.01,mmu:0.01):0.02,rno:0.03)"
    )
    writer = sql_writer(
        outpath, create=True, if_exists="overwrite" if overwrite else "raise"
    )
    app = renamer + aligner + writer
    app.apply_to(
        seq_collections,
        parallel=parallel,
        show_progress=verbose > 0,
        logger=LOGGER,
        cleanup=True,
    )
    result_dstore = io.get_data_store(outpath)
    rich_display(result_dstore.describe)
    if len(result_dstore.incomplete) > 0 and verbose:
        rich_display(result_dstore.summary_incomplete)

    return True


def filter_positions(inpath, outpath, overwrite, verbose):
    """filters alignment positions, removing non-canonical bases"""
    from cogent3.app import io, sample

    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    LOGGER.log_file_path = outpath.parent / "mdeqasis-fxy-filter_alignments.log"

    dstore = io.get_data_store(inpath)

    loader = sql_loader()

    just_nucs = sample.omit_degenerates(
        moltype="dna", motif_length=1, gap_is_degen=True
    )
    writer = sql_writer(
        outpath,
        create=True,
        if_exists="overwrite" if overwrite else "raise",
    )
    app = loader + just_nucs + writer
    r = app.apply_to(dstore, cleanup=True, show_progress=True, logger=LOGGER)
    app.data_store.close()
    result_dstore = io.get_data_store(outpath)
    rich_display(result_dstore.describe)
    if len(result_dstore.incomplete) > 0 and verbose:
        rich_display(result_dstore.summary_incomplete)

    return True
