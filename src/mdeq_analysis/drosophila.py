import re
from pathlib import Path

from cogent3 import get_app, open_data_store
from mdeq.utils import rich_display, summary_not_completed, write_to_sqldb
from mdeq_analysis.align_checked import alignment_filter, sequence_filter
from rich.progress import track
from scitrack import CachingLogger

_head = re.compile(r"##\s*(FBgn_ID|organism)")


def load_flybase_tsv(path, limit=None):
    from cogent3 import make_table, open_

    rows = []
    header = None
    num_cols = 0
    num_records = 0
    with open_(path) as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if header is None and not _head.search(line):
                continue

            if header is None:
                line = line.replace("##", "")

                header = [f.strip() for f in line.split("\t")]
                num_cols = len(header)
                continue

            if line.startswith("#"):  # comment line
                continue

            line = line.split("\t")
            if missing := num_cols - len(line):
                if missing > num_cols:
                    raise ValueError(f"too many columns {line}")

                line.extend([""] * missing)

            rows.append(line)
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
    dstore = open_data_store(indir, suffix="afa.mask", limit=limit)
    if testrun:
        print(dstore)
        return True

    reader = get_app("load_aligned", format="fasta", moltype="dna")
    take_seqs = get_app("take_named_seqs", "Dmel", "Dsim", "Dyak")
    align_qual = alignment_filter()
    seq_qual = sequence_filter()
    third_codons = get_app("take_codon_positions", 3)
    no_degenerates = get_app("omit_degenerates", moltype="dna", gap_is_degen=True)
    no_shorties = get_app("min_length", 300)

    if outpath.exists() and not overwrite:
        return outpath

    outpath.parent.mkdir(parents=True, exist_ok=True)
    out_dstore = open_data_store(outpath, mode="w")
    writer = write_to_sqldb(out_dstore)
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
    for m in track(dstore.completed, refresh_per_second=5):
        dmel_sym = m.unique_id.split("_")[0]
        if dmel_sym not in fb_ids:
            continue

        result = app(m)
        if not result:
            result.source = f"{dmel_sym}"
            if verbose > 2:
                print(result)
        else:
            result.info.source = f"{dmel_sym}"
            success += 1

        writer(result)

    log_file_path = Path(LOGGER.log_file_path)
    LOGGER.shutdown()
    out_dstore.write_log(unique_id=log_file_path.name, data=log_file_path.read_text())
    log_file_path.unlink()
    out_dstore.unlock()

    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    out_dstore.close()
    return True
