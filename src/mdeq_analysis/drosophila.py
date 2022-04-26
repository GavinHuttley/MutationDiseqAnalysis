import re

from pathlib import Path

from cogent3.app import io, sample
from mdeq.sqlite_data_store import sql_writer
from mdeq.utils import rich_display
from rich.progress import track
from scitrack import CachingLogger

from mdeq_analysis.align_checked import alignment_filter, sequence_filter


_head = re.compile("##\s*(FBgn_ID|organism)")


def load_flybase_tsv(path, limit=None):
    from cogent3 import make_table, open_

    rows = []
    header = None
    num_records = 0
    with open_(path) as infile:
        for line in infile:
            if not line.strip():
                continue
            if header is None and not _head.search(line):
                continue

            if header is None:
                line = line.replace("##", "")

                header = [f.strip() for f in line.split("\t")]
                continue
            if line.startswith("#"):  # comment line
                continue

            rows.append(line.split("\t"))
            num_records += 1
            if num_records == limit:
                break

    return make_table(header=header, data=rows, space=2)


def dros_sample_alignments(
    indir, fb_table, outpath, limit, overwrite, verbose, testrun
):
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    indir = Path(indir)
    outpath = Path(outpath)

    LOGGER.log_file_path = outpath.parent / "mdeqasis-dros_sample_alignments.log"
    LOGGER.input_file(fb_table)
    fb_table = load_flybase_tsv(fb_table)
    fb_ids = set(fb_table.columns["primary_FBid"])
    dstore = io.get_data_store(indir, suffix="afa.mask", limit=limit)
    if testrun:
        print(dstore)
        return True

    reader = io.load_aligned(format="fasta", moltype="dna")
    take_seqs = sample.take_named_seqs("Dmel", "Dsim", "Dyak")
    align_qual = alignment_filter()
    seq_qual = sequence_filter()
    third_codons = sample.take_codon_positions(3)
    no_degenerates = sample.omit_degenerates(moltype="dna", gap_is_degen=True)
    no_shorties = sample.min_length(300)
    writer = sql_writer(
        outpath,
        create=True,
        if_exists=io.OVERWRITE if overwrite else io.RAISE,
    )
    # we do the write as a separate step so that the source can be simplified
    app = (
        reader
        + take_seqs
        + align_qual
        + seq_qual
        + third_codons
        + no_degenerates
        + no_shorties
    )
    success = 0

    for m in track(dstore, refresh_per_second=5):
        dmel_sym = Path(m).stem.split("_")[0]
        if dmel_sym not in fb_ids:
            continue

        result = app(m)
        if not result:
            result.source = f"{dmel_sym}.json"
            if verbose > 2:
                print(result)
            writer.data_store.write_incomplete(result.source, result)
        else:
            result.info.source = f"{dmel_sym}.json"
            success += 1
            writer(result)

    log_file_path = LOGGER.log_file_path
    LOGGER.shutdown()
    writer.data_store.add_file(log_file_path, cleanup=True, keep_suffix=True)
    writer.data_store.close()

    rich_display(writer.data_store.describe)
    if len(writer.data_store.incomplete) > 0 and verbose:
        rich_display(writer.data_store.summary_incomplete)

    return True
