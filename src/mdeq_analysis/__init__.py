"""scripts for analysis steps associated with the mdeq app"""
import mdeq  # isort: skip  # make sure this stays at the top
import click


try:
    from wakepy import set_keepawake, unset_keepawake
except NotImplementedError:
    # probably on linux where this library doesn't work
    make_none = type(None)

    set_keepawake, unset_keepawake = make_none, make_none


from mdeq_analysis import microbial as micro


__author__ = "Gavin Huttley"
__credits__ = ["Gavin Huttley"]

__version__ = "2022.03.14"


@click.group()
def main():
    """cli for analyses assocuated with mdeq"""


_seed_alignment = click.option(
    "-sa",
    "--seed_alignment",
    type=click.Choice(["hi-hi", "hi-lo", "lo-hi", "lo-lo"]),
    required=True,
    help="specify seed alignment in terms of relative jsd, entropy",
)
_sim_length = click.option("-s", "--limit", type=int, default=None)
_indir = click.option(
    "-d",
    "--indir",
    type=click.Path(exists=True),
    help="directory containing input data",
)
_outdir = click.option(
    "-od",
    "--outdir",
    type=click.Path(),
    help="directory to write output",
)


@main.command()
@_indir
@click.option("-s", "--suffix", help="suffix of files within indir to be loaded")
@mdeq._outpath
@mdeq._limit
@mdeq._overwrite
def filter_alignments(**kwargs):
    """filters sequence alignments"""
    result = micro.filter_alignments(**kwargs)
    if result:
        click.secho("Done!", fg="green")
    else:
        click.secho("Failed!", fg="red")


@main.command()
@mdeq._inpath
@mdeq._outpath
@mdeq._parallel
@mdeq._mpi
@mdeq._limit
@mdeq._overwrite
@mdeq._verbose
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
@mdeq._inpath
@mdeq._outpath
@mdeq._parallel
@mdeq._limit
@mdeq._overwrite
@mdeq._verbose
def microbial_gn_stats(**kwargs):
    """generate stats from GN fits to microbial 16S data"""
    result = micro.gn_statistics(**kwargs)
    if result:
        click.secho("Done!", fg="green")
    else:
        click.secho("Failed!", fg="red")


if __name__ == "__main__":
    main()
