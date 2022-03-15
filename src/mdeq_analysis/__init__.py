"""scripts for analysis steps associated with the mdeq app"""
import mdeq  # isort: skip  # make sure this stays at the top
import click
from wakepy import set_keepawake, unset_keepawake

from mdeq_analysis import microbial as micro


__author__ = "Gavin Huttley"
__credits__ = ["Gavin Huttley"]

__version__ = "2022.03.14"


@click.group()
def main():
    """cli for analyses assocuated with mdeq"""


_inpath = click.option(
    "-i", "--inpath", type=click.Path(exists=True), help="path to a tinydb"
)
_outpath = click.option(
    "-o",
    "--outpath",
    type=click.Path(),
    help="path to create a result tinydb",
)

_verbose = click.option("-v", "--verbose", count=True)
_limit = click.option("-L", "--limit", type=int, default=None)
_overwrite = click.option("-O", "--overwrite", is_flag=True)
_parallel = click.option(
    "-p",
    "--parallel",
    is_flag=True,
    help="run in parallel (on single machine)",
)
_testrun = click.option(
    "-t",
    "--testrun",
    is_flag=True,
    help="don't write anything, quick (but inaccurate) optimisation",
)


@main.command()
@click.option(
    "-d",
    "--indir",
    type=click.Path(exists=True),
    help="directory containing input data",
)
@click.option("-s", "--suffix", help="suffix of files within indir to be loaded")
@click.option("-o", "--outpath", type=click.Path(), help="path to write tinydb output")
@_limit
@_overwrite
def filter_alignments(**kwargs):
    """filters sequence alignments"""
    result = micro.filter_alignments(**kwargs)
    if result:
        click.secho("Done!", fg="green")
    else:
        click.secho("Failed!", fg="red")


@main.command()
@_inpath
@_outpath
@_parallel
@_limit
@_overwrite
@_verbose
def microbial_fit_gn(**kwargs):
    """fits GN to microbial 16S data"""
    set_keepawake(keep_screen_awake=False)
    result = micro.fit_gn(**kwargs)
    if result:
        click.secho("Done!", fg="green")
    else:
        click.secho("Failed!", fg="red")
    unset_keepawake()


@main.command()
@_inpath
@_outpath
@_parallel
@_limit
@_overwrite
@_verbose
def microbial_gn_stats(**kwargs):
    """generate stats from GN fits to microbial 16S data"""
    result = micro.gn_statistics(**kwargs)
    if result:
        click.secho("Done!", fg="green")
    else:
        click.secho("Failed!", fg="red")


if __name__ == "__main__":
    main()
