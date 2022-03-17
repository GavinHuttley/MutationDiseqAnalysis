"""sampling and analysis of microbial data"""
import pathlib

from dataclasses import asdict, dataclass
from random import Random

from cogent3.app.composable import SERIALISABLE_TYPE, appify, get_data_source
from cogent3.util import deserialise, misc
from mdeq.jsd import get_entropy, get_jsd
from mdeq.model import mles_within_bounds
from mdeq.stationary_pi import get_stat_pi_via_eigen
from mdeq.utils import (
    SerialisableMixin,
    configure_parallel,
    foreground_from_jsd,
)
from numpy.linalg import cond, eig
from scitrack import CachingLogger


__author__ = "Gavin Huttley"
__credits__ = ["Kath Caley", "Gavin Huttley"]

__version__ = "2022.03.14"


def _make_name(x):
    p = get_data_source(x)
    return pathlib.Path(pathlib.Path(p).stem).stem


def filter_alignments(indir, outpath, suffix, limit, overwrite):
    """filters sequence alignments"""
    from cogent3.app import io, sample

    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    indir = pathlib.Path(indir)
    outpath = pathlib.Path(outpath)

    LOGGER.log_file_path = outpath.parent / "mdeqasis-filter_alignments.log"

    dstore = io.get_data_store(indir, suffix=suffix, limit=limit)

    loader = io.load_aligned(moltype="dna", format="nexus")

    just_nucs = sample.omit_degenerates(
        moltype="dna", motif_length=1, gap_is_degen=True
    )
    writer = io.write_db(
        outpath,
        name_callback=_make_name,
        create=True,
        if_exists="overwrite" if overwrite else "raise",
    )
    app = loader + just_nucs + writer
    r = app.apply_to(dstore, cleanup=True, show_progress=True, logger=LOGGER)
    app.data_store.close()
    return True


def fit_gn(inpath, outpath, parallel, mpi, limit, overwrite, verbose):
    """fits GN to microbial 16S data"""
    from cogent3.app import evo, io

    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    inpath = pathlib.Path(inpath)
    outpath = pathlib.Path(outpath)

    LOGGER.log_file_path = outpath.parent / "mdeqasis-fit_gn.log"

    dstore = io.get_data_store(inpath, limit=limit)
    overwrite = "overwrite" if overwrite else "raise"
    loader = io.load_db()
    gn = evo.model(
        "GN",
        unique_trees=True,
        time_het="max",
        optimise_motif_probs=True,
        opt_args=dict(max_restarts=5),
        show_progress=verbose > 2,
    )
    writer = io.write_db(outpath, create=True, if_exists=overwrite)

    app = loader + gn + writer
    dstore = io.get_data_store(inpath, limit=limit)

    kwargs = configure_parallel(parallel=parallel, mpi=mpi)

    r = app.apply_to(
        dstore, cleanup=True, logger=LOGGER, show_progress=verbose >= 2, **kwargs
    )

    print(app.data_store.describe)
    app.data_store.close()
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
    def header(self):
        return self.__slots__


@deserialise.register_deserialiser(misc.get_object_provenance(gn_fit_stats))
def deserialise_gn_fit_stats(data):
    return gn_fit_stats.from_dict(data)


@appify(input_types=SERIALISABLE_TYPE, output_types=SERIALISABLE_TYPE)
def compute_stats(result):
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
    edge = foreground_from_jsd(aln)
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
    from cogent3.app import io

    LOGGER = CachingLogger(create_dir=True)
    LOGGER.log_args()

    inpath = pathlib.Path(inpath)
    outpath = pathlib.Path(outpath)

    LOGGER.log_file_path = outpath.parent / "mdeqasis-gn_statistics.log"

    dstore = io.get_data_store(inpath, limit=limit)

    overwrite = "overwrite" if overwrite else "raise"
    loader = io.load_db()
    calc_stats = compute_stats()
    writer = io.write_db(outpath, create=True, if_exists=overwrite)
    within_bounds = mles_within_bounds()
    app = loader + within_bounds + calc_stats + writer
    r = app.apply_to(
        dstore,
        logger=LOGGER,
        cleanup=True,
        show_progress=verbose > 1,
        parallel=parallel,
    )
    print(app.data_store.describe)
    app.data_store.close()
    return True
