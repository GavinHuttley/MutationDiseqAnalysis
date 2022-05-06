"""scripts for analysis steps associated with the mdeq app"""
import mdeq  # isort: skip  # make sure this stays at the top
import glob
import inspect

from pathlib import Path
from warnings import simplefilter

import click
import mdeq._click_options as mdeq_opt


try:
    from wakepy import set_keepawake, unset_keepawake
except NotImplementedError:
    # probably on linux where this library doesn't work
    make_none = type(None)

    set_keepawake, unset_keepawake = make_none, make_none


from mdeq_analysis import fxy
from mdeq_analysis import microbial as micro
from mdeq_analysis import postprocess


simplefilter("ignore", category=UserWarning, append=True)

__author__ = "Gavin Huttley"
__credits__ = ["Gavin Huttley"]

__version__ = "2022.03.14"


@click.group()
def main():
    """cli for analyses associated with mdeq"""


_seed_aln = click.option(
    "-sa",
    "--seed_aln",
    type=click.Choice(["hi_hi", "hi_lo", "lo_hi", "lo_lo"]),
    required=True,
    help="specify seed alignment in terms of relative jsd, entropy",
)
_sim_length = click.option(
    "-sl",
    "--sim_length",
    type=click.Choice(["300", "3000", "30000"]),
    required=True,
    help="alignment length to generate",
)
_indir = click.option(
    "-d",
    "--indir",
    type=click.Path(exists=True),
    help="directory containing input data",
)


@main.command()
@_indir
@click.option("-s", "--suffix", help="suffix of files within indir to be loaded")
@mdeq_opt._outpath
@mdeq_opt._limit
@mdeq_opt._overwrite
def filter_alignments(**kwargs):
    """filters sequence alignments"""
    result = micro.filter_alignments(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


@main.command()
@mdeq_opt._inpath
@mdeq_opt._outpath
@mdeq_opt._parallel
@mdeq_opt._mpi
@mdeq_opt._limit
@mdeq_opt._overwrite
@mdeq_opt._verbose
def microbial_fit_gn(**kwargs):
    """fits GN to microbial 16S data"""
    set_keepawake(keep_screen_awake=False)
    result = micro.fit_gn(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")
    unset_keepawake()


@main.command()
@mdeq_opt._inpath
@mdeq_opt._outpath
@mdeq_opt._parallel
@mdeq_opt._limit
@mdeq_opt._overwrite
@mdeq_opt._verbose
def microbial_gn_stats(**kwargs):
    """generate stats from GN fits to microbial 16S data"""
    result = micro.gn_statistics(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


@main.command()
@mdeq_opt._inpath
@mdeq_opt._outdir
@mdeq_opt._just_continuous
@_seed_aln
@mdeq_opt._seed
@_sim_length
@mdeq_opt._num_reps
@mdeq_opt._overwrite
@mdeq_opt._verbose
@mdeq_opt._testrun
def microbial_fg_gsn_synthetic(**kwargs):
    """generate synthetic alignments from GN fits to microbial 16S data"""
    result = micro.GSN_synthetic(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


def _all_tinydbs(*args):
    indir = args[-1]
    if "*" not in indir:
        return list(Path(indir).glob("*.tinydb"))

    return [Path(p) for p in glob.glob(indir) if p.endswith("tinydb")]


_glob_indir = click.option(
    "-d",
    "--indir",
    callback=_all_tinydbs,
    help="more general path allowing glob patterns",
)


@main.command()
@_glob_indir
@mdeq_opt._outdir
@mdeq_opt._wrt_nstat
@mdeq_opt._limit
@mdeq_opt._overwrite
@mdeq_opt._verbose
def microbial_nabla(**kwargs):
    """generate nabla statistics from GN fits of TOE boostrap results on null data"""
    paths = kwargs.pop("indir")
    func_name = inspect.stack()[0].function
    if kwargs["verbose"]:
        print(paths)

    set_keepawake(keep_screen_awake=False)

    for i, path in enumerate(
        map(lambda x: micro.generate_convergence(x, **kwargs), paths)
    ):
        if path:
            click.secho(f"{func_name!r} success for {paths[i].name!r}", fg="green")
        else:
            click.secho(f"{func_name!r} failed  for {paths[i].name!r}", fg="red")

    unset_keepawake()


@main.command()
@_glob_indir
@mdeq_opt._limit
@mdeq_opt._overwrite
@mdeq_opt._verbose
@mdeq_opt._testrun
def microbial_aeop_locations(**kwargs):
    """generates a synthetic locations file for alignments at inpath"""
    func_name = inspect.stack()[0].function
    paths = kwargs.pop("indir")
    for path in paths:
        result = micro.make_synthetic_aeop_locations(path, **kwargs)
        if result:
            click.secho(f"{func_name!r} wrote {result!r}!", fg="green")
        else:
            click.secho(f"{func_name!r} failed for {path!r}!", fg="red")

        if kwargs["testrun"]:
            break


@main.command()
@mdeq_opt._inpath
@mdeq_opt._outdir
@_seed_aln
@mdeq_opt._seed
@_sim_length
@mdeq_opt._num_reps
@mdeq_opt._overwrite
@mdeq_opt._verbose
@mdeq_opt._testrun
def microbial_teop_synthetic(**kwargs):
    """generate synthetic alignments from GN fits to microbial 16S data"""
    result = micro.make_synthetic_teop(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


@main.command()
@_indir
@click.option(
    "-fb", "--fb_table", required=True, type=Path, help="flybase gene order table"
)
@mdeq_opt._outpath
@mdeq_opt._limit
@mdeq_opt._overwrite
@mdeq_opt._verbose
@mdeq_opt._testrun
def dros_filter_alignments(**kwargs):
    """filters Drosophila group sequence alignments"""
    from mdeq_analysis import drosophila as dros

    result = dros.dros_sample_alignments(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


@main.command()
@click.option(
    "-cd",
    "--cds_indir",
    required=True,
    type=Path,
    help="directory of unaligned CDS sequences",
)
@click.option(
    "-id",
    "--intron_indir",
    required=True,
    type=Path,
    help="directory of ensembl aligned intron sequences",
)
@click.option(
    "-md", "--metadata_path", required=True, type=Path, help="tsv gene order table"
)
@mdeq_opt._outpath
@mdeq_opt._limit
@mdeq_opt._overwrite
@mdeq_opt._verbose
@mdeq_opt._testrun
def ape_align_cds(**kwargs):
    """codon aligns chrom1 human, chimp, gorilla sequences"""
    from mdeq_analysis import ape

    result = ape.codon_align(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


@main.command()
@click.option(
    "-cp",
    "--cds_path",
    required=True,
    callback=mdeq_opt._valid_sqlitedb_input,
    help="directory of unaligned CDS sequences",
)
@click.option(
    "-id",
    "--intron_indir",
    required=True,
    type=Path,
    help="directory of ensembl aligned intron sequences",
)
@mdeq_opt._outdir
@mdeq_opt._limit
@mdeq_opt._overwrite
@mdeq_opt._verbose
@mdeq_opt._testrun
def ape_match_cds_intron(**kwargs):
    """filters ape codon and intron alignments, only writes out if gene passes both filters"""
    from mdeq_analysis import ape

    result = ape.match_cds_intron(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


@main.command()
@click.option(
    "-mmu",
    "--mmu_path",
    required=True,
    type=Path,
    help="directory of unaligned Mus musculus sequences",
)
@click.option(
    "-msp",
    "--msp_path",
    required=True,
    type=Path,
    help="directory of unaligned Mus spretus sequences",
)
@click.option(
    "-rno",
    "--rno_path",
    required=True,
    type=Path,
    help="directory of unaligned Rattus norvegicus sequences",
)
@mdeq_opt._outpath
@mdeq_opt._overwrite
@mdeq_opt._verbose
@mdeq_opt._parallel
def fxy_align(*args, **kwargs):
    """align rodent sequences"""
    result = fxy.make_aligned(*args, **kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


@main.command()
@mdeq_opt._inpath
@mdeq_opt._outpath
@mdeq_opt._overwrite
@mdeq_opt._verbose
def fxy_filter_positions(**kwargs):
    """filter Fxy sequence alignment positions, removing non-canonical bases"""
    result = fxy.filter_positions(**kwargs)
    func_name = inspect.stack()[0].function
    if result:
        click.secho(f"{func_name!r} is done!", fg="green")
    else:
        click.secho(f"{func_name!r} failed!", fg="red")


if __name__ == "__main__":
    main()
