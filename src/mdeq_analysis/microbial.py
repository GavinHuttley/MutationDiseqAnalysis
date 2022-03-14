"""sampling and analysis of microbial data"""
import pathlib

from dataclasses import dataclass

from cogent3.app.composable import SERIALISABLE_TYPE, appify
from cogent3.util import deserialise, misc
from mdeq.jsd import get_entropy, get_jsd
from mdeq.utils import SerialisableMixin, foreground_from_jsd
from numpy.linalg import cond, eig
from scitrack import CachingLogger


__author__ = "Gavin Huttley"
__credits__ = ["Kath Caley", "Gavin Huttley"]

__version__ = "2022.03.14"


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
        outpath, create=True, if_exists="overwrite" if overwrite else "raise"
    )
    app = loader + just_nucs + writer
    r = app.apply_to(dstore, cleanup=True, show_progress=True, logger=LOGGER)
    return True


def fit_gn(inpath, outpath, parallel, limit, overwrite, verbose):
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
        optimise_motif_probs=True,
        opt_args=dict(max_restarts=5),
    )
    writer = io.write_db(outpath, create=True, if_exists=overwrite)

    app = loader + gn + writer
    dstore = io.get_data_store(inpath, limit=limit)
    r = app.apply_to(
        dstore,
        cleanup=True,
        parallel=parallel,
        logger=LOGGER,
        show_progress=verbose >= 2,
    )
    print("", app.data_store.describe, "", sep="\n")
    app.data_store.close()
    return True


@dataclass
class gn_fit_stats(SerialisableMixin):
    source: str
    foreground: str
    jsd: float
    entropy: float
    cond_num: float


@deserialise.register_deserialiser(misc.get_object_provenance(gn_fit_stats))
def deserialise_gn_fit_stats(data):
    return gn_fit_stats.from_dict(data)


@appify(input_types=SERIALISABLE_TYPE, output_types=SERIALISABLE_TYPE)
def compute_stats(model_result):
    """
    Parameters
    ----------
    model_result : cogent3.app.result.model_result

    Returns
    -------
    JSD, fg_edge, entropy, and condition number of a model_result as a generic_result,
    for wrapping in a user_function.

    """
    aln = model_result.alignment
    edge = foreground_from_jsd(aln)
    edge, ingroup, jsd = get_jsd(aln)

    entropy = get_entropy(model_result, edge, stat_pi=True)

    Q = model_result.lf.get_rate_matrix_for_edge(edge)
    ev = eig(Q)[1]

    cond_num = cond(ev)

    return gn_fit_stats(
        source=str(aln.info.source),
        jsd=jsd,
        entropy=entropy,
        cond_num=cond_num,
        foreground=ingroup,
    )


def gn_statistics(inpath, outpath, limit, overwrite, verbose):
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
    app = loader + calc_stats + writer
    r = app.apply_to(dstore, logger=LOGGER, cleanup=True, show_progress=verbose > 1)
    print(app.data_store.describe)
    return True
