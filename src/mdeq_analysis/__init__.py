"""scripts for analysis steps associated with the mdeq app"""

import click
from scitrack import CachingLogger

__author__ = "Gavin Huttley"
__credits__ = ["Gavin Huttley"]

__version__ = "2022.03.14"


@click.group()
def main():
    """cli for analyses assocuated with mdeq"""


@main.command()
@click.option(
    "-d",
    "--indir",
    type=click.Path(exists=True),
    help="directory containing input data",
)
@click.option("-s", "--suffix", help="suffix of files within indir to be loaded")
@click.option("-o", "--outpath", type=click.Path(), help="path to write tinydb output")
def filter_alignments(indir, outpath, suffix):
    """filters sequence alignments"""
    from cogent3.app import io, sample

    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()
    LOGGER.log_file_path = outpath.parent / "mdeq_analysis-filter.log"

    path = "../data/raw/microbial/nexus"
    dstore = io.get_data_store(path, suffix=suffix)

    loader = io.load_aligned(moltype="dna", format="nexus")

    just_nucs = sample.omit_degenerates(
        moltype="dna", motif_length=1, gap_is_degen=True
    )
    outpath = "../data/raw/microbial/filtered.tinydb"
    writer = io.write_db(outpath, create=True, if_exists="overwrite")
    app = loader + just_nucs + writer
    r = app.apply_to(dstore, cleanup=True, show_progress=True, logger=LOGGER)
    click.secho("Done!", fg="green")


if __name__ == "__main__":
    main()
