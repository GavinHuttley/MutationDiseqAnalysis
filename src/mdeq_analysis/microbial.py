"""sampling and analysis of microbial data"""
import inspect
import pathlib

from dataclasses import asdict, dataclass

from cogent3 import get_app, make_table, open_data_store
from cogent3.app import typing as c3_types
from cogent3.app.composable import define_app, get_data_source
from cogent3.util import deserialise, misc
from mdeq.jsd import get_entropy, get_jsd
from mdeq.model import mles_within_bounds
from mdeq.stationary_pi import get_stat_pi_via_eigen
from mdeq.utils import (
    SerialisableMixin,
    configure_parallel,
    foreground_from_jsd,
    load_from_sqldb,
    rich_display,
    summary_not_completed,
    write_to_sqldb,
)
from numpy import iinfo, int64, random
from numpy.linalg import cond, eig
from rich.progress import track
from scitrack import CachingLogger


__author__ = "Gavin Huttley"
__credits__ = ["Kath Caley", "Gavin Huttley"]
__version__ = "2022.03.14"


def load_seed_alignment(inpath, seed_aln):
    """loads the alignment with matching name to seed_aln"""
    dstore = open_data_store(inpath)
    r = [m for m in dstore if seed_aln in m.unique_id]
    if len(r) == 0:
        raise ValueError(f"{seed_aln=!r} missing from {inpath!r}")
    loader = load_from_sqldb()
    return loader(r[0])


def _make_name(x):
    p = get_data_source(x)
    return pathlib.Path(pathlib.Path(p).stem).stem


def filter_alignments(indir, outpath, suffix, limit, overwrite):
    """filters sequence alignments"""
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    indir = pathlib.Path(indir)
    outpath = pathlib.Path(outpath)

    LOGGER.log_file_path = outpath.parent / "mdeqasis-filter_alignments.log"

    dstore = open_data_store(indir, suffix=suffix, limit=limit)

    loader = get_app("load_aligned", moltype="dna", format="nexus")

    just_nucs = get_app(
        "omit_degenerates", moltype="dna", motif_length=1, gap_is_degen=True
    )

    if outpath.exists() and not overwrite:
        raise IOError(f"{str(outpath)} exists")

    out_dstore = open_data_store(outpath, mode="w")

    writer = write_to_sqldb(
        out_dstore,
        id_from_source=_make_name,
    )
    app = loader + just_nucs + writer
    out_dstore = app.apply_to(
        dstore.completed, cleanup=True, show_progress=True, logger=LOGGER
    )
    out_dstore.unlock()

    rich_display(app.data_store.describe)
    if len(app.data_store.not_completed) > 0:
        rich_display(summary_not_completed(out_dstore))

    out_dstore.close()

    return True


def fit_gn(inpath, outpath, parallel, mpi, limit, overwrite, verbose):
    """fits GN to microbial 16S data"""
    from cogent3.app import evo, io

    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    inpath = pathlib.Path(inpath)
    outpath = pathlib.Path(outpath)

    LOGGER.log_file_path = outpath.parent / "mdeqasis-fit_gn.log"

    loader = load_from_sqldb()
    gn = get_app(
        "model",
        "GN",
        unique_trees=True,
        time_het="max",
        optimise_motif_probs=True,
        opt_args=dict(max_restarts=5),
        show_progress=verbose > 2,
    )

    if outpath.exists() and not overwrite:
        raise IOError(f"{str(outpath)} exists")

    out_dstore = open_data_store(outpath, mode="w")

    writer = write_to_sqldb(out_dstore)

    app = loader + gn + writer
    dstore = open_data_store(inpath, limit=limit)

    kwargs = configure_parallel(parallel=parallel, mpi=mpi)

    r = app.apply_to(
        dstore.completed,
        cleanup=True,
        logger=LOGGER,
        show_progress=verbose >= 2,
        **kwargs,
    )

    out_dstore.unlock()

    rich_display(app.data_store.describe)
    if len(app.data_store.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    app.data_store.close()

    return True


def select_alignments(aligns_path, fit_path, outpath, limit, overwrite, verbose):
    """selects alignments based on GN fits"""
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()
    LOGGER.log_file_path = outpath.parent / "mdeqasis-micro_select_alignments.log"

    dstore = open_data_store(fit_path)
    loader = load_from_sqldb()
    records = [loader(m) for m in dstore.completed]
    header = records[0].header()
    rows = [r.to_record() for r in records]
    table = make_table(header=header, data=rows)
    table = table.filtered(lambda x: x <= 2, columns="cond_num")
    source = [f"{s}.json" for s in table.columns["source"]]
    align_dstore = open_data_store(aligns_path)
    sampled = [m for m in align_dstore.completed if m.unique_id in source]

    if outpath.exists() and not overwrite:
        raise IOError(f"{str(outpath)} exists")

    out_dstore = open_data_store(outpath, mode="w")

    writer = write_to_sqldb(out_dstore)
    for m in track(sampled):
        aln = loader(m)
        aln.info.source = m.unique_id
        writer(aln)

    log_file_path = pathlib.Path(LOGGER.log_file_path)
    LOGGER.shutdown()

    out_dstore.write_log(unique_id=log_file_path.name, data=log_file_path.read_text())
    log_file_path.unlink(missing_ok=True)
    out_dstore.unlock()

    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    out_dstore.close()
    return True


@dataclass(slots=True)
class gn_fit_stats(SerialisableMixin):
    source: str
    foreground: str
    jsd: float
    entropy: float
    cond_num: float

    def to_dict(self) -> dict:
        return asdict(self)

    def to_record(self):
        return tuple(self.to_dict().values())

    @classmethod
    def header(cls):
        return cls.__slots__


@deserialise.register_deserialiser(misc.get_object_provenance(gn_fit_stats))
def deserialise_gn_fit_stats(data):
    return gn_fit_stats.from_dict(data)


@define_app
def compute_stats(result: c3_types.SerialisableType) -> c3_types.SerialisableType:
    """
    Parameters
    ----------
    model_result : cogent3.app.result.model_result

    Returns
    -------
    JSD, fg_edge, entropy, and condition number of a model_result as a generic_result,
    for wrapping in a user_function.

    """
    if not result:  # handle case where result is NotCompleted
        return result

    aln = result.alignment
    edge, _, jsd = get_jsd(aln)

    entropy = get_entropy(result, edge, stat_pi=True)

    Q = result.lf.get_rate_matrix_for_edge(edge)
    ev = eig(Q)[1]

    cond_num = cond(ev)

    return gn_fit_stats(
        source=_make_name(aln.info.source),
        jsd=jsd,
        entropy=float(entropy),
        cond_num=cond_num,
        foreground=edge,
    )


def gn_statistics(inpath, outpath, parallel, limit, overwrite, verbose):
    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    inpath = pathlib.Path(inpath)
    outpath = pathlib.Path(outpath)

    LOGGER.log_file_path = outpath.parent / "mdeqasis-gn_statistics.log"

    dstore = open_data_store(inpath, limit=limit)

    loader = load_from_sqldb()
    calc_stats = compute_stats()

    if outpath.exists() and not overwrite:
        raise IOError(f"{str(outpath)} exists")

    outpath.parent.mkdir(exist_ok=True, parents=True)
    out_dstore = open_data_store(outpath, mode="w")

    writer = write_to_sqldb(out_dstore)

    within_bounds = mles_within_bounds()
    app = loader + within_bounds + calc_stats + writer
    out_dstore = app.apply_to(
        dstore.completed,
        logger=LOGGER,
        cleanup=True,
        show_progress=verbose > 1,
        parallel=parallel,
    )

    out_dstore.unlock()

    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    out_dstore.close()

    return True


# the selected seed alignments, keys are `jsd-entropy`

seed_alignments = [
    ("hi_hi", "197113_332182_17210"),
    ("lo_hi", "198257_206396_13724"),
    ("hi_lo", "200580_114946_573911"),
    ("lo_lo", "758_443154_73021"),
]


def GSN_synthetic(
    inpath,
    outdir,
    just_continuous,
    seed_aln,
    seed,
    sim_length,
    num_reps,
    overwrite,
    verbose,
    testrun,
):
    """simulate alignments under mixed general and stationary models

    Parameters
    ----------
    inpath : str
        data store containing observed alignment
    outdir : str
        dir to write output
    just_continuous : bool
        if True, only continuous-time process across entire tree and a GSN
        model will be fit
    seed_aln : str
        name of alignment to be used as seed
    seed : int
        random number seed
    sim_length : int
        length of simulated alignment
    num_reps : int
        number of alignments to simulate
    """
    LOGGER = CachingLogger(create_dir=True)

    inpath = pathlib.Path(inpath)
    outdir = pathlib.Path(outdir)
    func_name = inspect.stack()[0].function
    if just_continuous:
        outpath = (
            outdir / f"{func_name}-{seed_aln}-{sim_length}bp-{num_reps}reps.sqlitedb"
        )
    else:
        outpath = (
            outdir / f"fg_{func_name}-{seed_aln}-{sim_length}bp-{num_reps}reps.sqlitedb"
        )

    outpath.parent.mkdir(exist_ok=True, parents=True)

    LOGGER.log_args()

    LOGGER.log_file_path = f"{outpath.stem}.log"
    LOGGER.input_file(inpath)

    obs_aln = load_seed_alignment(inpath, dict(seed_alignments)[seed_aln])
    # model foreground means has discrete edges for background
    has_discrete = not just_continuous
    fg_edge, lf = _get_null_generator(obs_aln, has_discrete, testrun, verbose)

    if verbose > 1:
        print(lf)

    if outpath.exists() and not overwrite:
        raise IOError(f"{str(outpath)} exists")

    out_dstore = open_data_store(outpath, mode="w")

    writer = write_to_sqldb(out_dstore)

    # make a seeded rng
    rng = random.default_rng(seed=seed)
    # we will use numpy to select a new seed each iteration, so we identify
    # the upper limit to choose from as the maximum for a 64-bit integer
    max_int = iinfo(int64).max
    sim_seeds = set()
    sim_length = int(sim_length)
    for i in track(range(num_reps)):
        while True:
            # ensure the seed is unique
            sim_seed = rng.choice(max_int)
            if sim_seed not in sim_seeds:
                sim_seeds.add(sim_seed)
                break

        sim_aln = lf.simulate_alignment(sequence_length=sim_length, seed=sim_seed)
        sim_aln.info.source = f"{seed_aln}-sim-{i}.json"
        if has_discrete:
            sim_aln.info.fg_edge = fg_edge

        writer(sim_aln)

    log_file_path = pathlib.Path(LOGGER.log_file_path)
    LOGGER.shutdown()
    out_dstore.write_log(unique_id=log_file_path.name, data=log_file_path.read_text())
    log_file_path.unlink()

    out_dstore.unlock()

    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    out_dstore.close()

    return True


def _get_null_generator(obs_aln, has_discrete, testrun, verbose):
    # figure out the fg edge
    if has_discrete:
        model_name = "GN"
        fg_edge, _, _ = get_jsd(obs_aln)
        bg_edges = list({fg_edge} ^ set(obs_aln.names))
        lf_args = dict(discrete_edges=bg_edges, expm="pade")
    else:
        model_name = "GSN"
        lf_args = dict(expm="pade")
        fg_edge = None

    if testrun:
        opt_args = {"max_restarts": 1, "max_evaluations": 10, "limit_action": "ignore"}
    else:
        opt_args = {"max_restarts": 5, "tolerance": 1e-8}

    gn = get_app(
        "model",
        model_name,
        optimise_motif_probs=True,
        lf_args=lf_args,
        opt_args=opt_args,
        show_progress=verbose > 2,
    )
    lf = gn(obs_aln).lf
    if not has_discrete:
        return None, lf

    psub_fg = lf.get_psub_for_edge(fg_edge)
    stat_pi_fg = get_stat_pi_via_eigen(psub_fg)
    # these become the root probabilities so the fg edge is stationary in
    # the synthetic data
    pi = {base: stat_pi_fg[index] for index, base in enumerate(obs_aln.moltype)}
    lf.set_motif_probs(pi)
    return fg_edge, lf


def generate_convergence(
    path, outdir, wrt_nstat, limit=None, overwrite=False, verbose=0
):
    """produces nabla statistics from TOE results"""
    from mdeq.convergence import bootstrap_to_nabla

    name_suffix = (
        "nstat_scale-convergence.sqlitedb"
        if wrt_nstat
        else "standard_scale-convergence.sqlitedb"
    )
    outpath = outdir / f"{path.stem}-{name_suffix}"
    if outpath.exists() and not overwrite:
        return outpath

    dstore = open_data_store(path, limit=limit)
    limit = limit or 1000
    if len(dstore) != limit:
        print(f"{path.stem} has {limit - len(dstore)} incomplete results!")
        summary_not_completed(dstore)
        if len(dstore.completed) == 0:
            return False

    if verbose:
        print("getting nabla values")

    loader = load_from_sqldb()
    conv = bootstrap_to_nabla()

    outpath.parent.mkdir(exist_ok=True, parents=True)
    out_dstore = open_data_store(outpath, mode="w")
    writer = write_to_sqldb(out_dstore)

    app = loader + conv + writer
    out_dstore = app.apply_to(
        dstore.completed,
        cleanup=True,
        show_progress=verbose > 1,
    )

    out_dstore.unlock()
    success = len(out_dstore.completed) > 0
    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    app.data_store.close()

    return success


def make_synthetic_aeop_locations(
    inpath,
    limit,
    overwrite,
    verbose,
    testrun,
):
    """creates synthetic locations suitable for input into mdeq aeop
    for the alignments at inpath

    Notes
    -----
    To ensure independence of analysed results, the locations file
    assigns 2 alignments to a unique 'chromosome'!
    Replaces the sqlitedb suffix for tsv.
    """
    outpath = inpath.parent / pathlib.Path(f"locations-{inpath.stem}.tsv")
    if outpath.exists() and not overwrite:
        return outpath

    dstore = open_data_store(inpath, limit=limit)
    loader = load_from_sqldb()
    names = [pathlib.Path(loader(m).info.source).stem for m in dstore]
    if verbose:
        print(names[:6])

    records = []
    for i in range(0, len(names), 2):
        pair = names[i : i + 2]
        records.extend([[n, i + 1, j] for j, n in enumerate(pair)])
    table = make_table(header=["name", "coord_name", "start"], data=records)
    if verbose:
        print(table[:6])

    if not testrun:
        table.write(outpath)

    return outpath


seed_outgroup = {
    "hi_hi": "197113",
    "hi_lo": "200580",
    "lo_hi": "198257",
    "lo_lo": "758",
}


def make_synthetic_teop(
    inpath,
    outdir,
    seed_aln,
    seed,
    sim_length,
    num_reps,
    overwrite,
    verbose,
    testrun,
):
    """creates synthetic alignments under teop null from seed alignment

    Notes
    -----
    The null model is one in which the ingroup have the same process, which
    differs from the outgroup.
    """
    from mdeq.eop import temporal_eop

    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    outpath = outdir / pathlib.Path(
        f"teop-{seed_aln}-{sim_length}bp-{num_reps}reps.sqlitedb"
    )
    if outpath.exists() and not overwrite:
        return outpath

    LOGGER.log_file_path = f"{outpath.stem}.log"
    LOGGER.input_file(inpath)

    seed_name = dict(seed_alignments)[seed_aln]
    obs_aln = load_seed_alignment(inpath, seed_name)

    ingroup = set(obs_aln.names) ^ {seed_outgroup[seed_aln]}

    teop = temporal_eop(
        ingroup,
    )
    if outpath.exists() and not overwrite:
        raise IOError(f"{str(outpath)} exists")

    outpath.parent.mkdir(parents=True, exist_ok=True)
    out_dstore = open_data_store(outpath, mode="w")

    writer = write_to_sqldb(out_dstore)

    # make a seeded rng
    rng = random.default_rng(seed=seed)
    # we will use numpy to select a new seed each iteration, so we identify
    # the upper limit to choose from as the maximum for a 64-bit integer
    max_int = iinfo(int64).max
    sim_seeds = set()

    sim_length = int(sim_length)
    fitted = teop(obs_aln)
    simulator = fitted.null.lf.simulate_alignment
    for i in track(range(num_reps)):
        while True:
            # ensure the seed is unique
            sim_seed = rng.choice(max_int)
            if sim_seed not in sim_seeds:
                sim_seeds.add(sim_seed)
                break

        sim_aln = simulator(sequence_length=sim_length, seed=sim_seed)
        sim_aln.info.source = f"{seed_name}-sim-{i}.json"
        writer(sim_aln)

    log_file_path = pathlib.Path(LOGGER.log_file_path)
    LOGGER.shutdown()
    out_dstore.write_log(unique_id=log_file_path.name, data=log_file_path.read_text())
    log_file_path.unlink()

    out_dstore.unlock()

    rich_display(out_dstore.describe)
    if len(out_dstore.not_completed) > 0 and verbose:
        rich_display(summary_not_completed(out_dstore))

    out_dstore.close()

    return True
