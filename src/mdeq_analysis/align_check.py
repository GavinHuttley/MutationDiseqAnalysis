import numpy as np

from cogent3.app.composable import (
    ALIGNED_TYPE,
    SEQUENCE_TYPE,
    SERIALISABLE_TYPE,
    NotCompleted,
    appify,
)


__author__ = ["Gavin Huttley", "Chris Bradley"]
__credits__ = ["Chris Bradley", "Gavin Huttley", "Kath Caley"]
__version__ = "2022.03.14"

_types = ALIGNED_TYPE, SERIALISABLE_TYPE


@appify(_types, _types)
def alignment_filter(align, threshold=1.5):
    """filter out poorly aligned sequences"""
    entropy = align.entropy_per_pos()
    entropy = np.where(np.isnan(entropy), 0, entropy)
    n_groups = len(entropy) / 21
    entropy_groups = np.array_split(entropy, n_groups)
    entropy_means = [np.mean(a) for a in entropy_groups]
    std = np.std(entropy_means)
    mean = np.mean(entropy_means)
    cev = std / mean
    if std != 0 and cev >= threshold:
        return NotCompleted(
            "FAIL", "alignment_filter", f"threshold < cev {cev} ", align
        )

    align.info["alignment_cev"] = cev
    return align


@appify(_types, _types)
def sequence_filter(align, threshold=50):
    """filters an alignment with a sequence with lots of repeats

    Parameters
    ----------
    align : Alignment
         sequence alignment
    threshold : int, optional
        coefficient of variation threshold, any sequence
        exceeding this (default 50) causes the entire alignment to
        to be excluded

    Notes
    -----
    Examines entropies for k in {3, 6, 9, 12} and determines the
    coefficient of variation (cev) for each. The mean and standard deviation
    of these entropies are used to compute the cev.
    """
    if not align:
        return align

    ks = [3, 6, 9, 12]
    entropies = [align.entropy_per_seq(k) for k in ks]
    std = np.std(entropies)
    mean = np.mean(entropies)
    cev = std / mean
    if std != 0 and cev >= threshold:
        return NotCompleted(
            "FAIL",
            "sequence_filter",
            f"threshold < cev {cev} ",
            source=align.info.source,
        )

    align.info["sequence_cev"] = cev
    return align
