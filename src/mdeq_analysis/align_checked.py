import numpy as np

from cogent3.app.composable import define_app, NotCompleted
from cogent3.app import typing as c3_types


__author__ = ["Gavin Huttley", "Chris Bradley"]
__credits__ = ["Chris Bradley", "Gavin Huttley", "Kath Caley"]
__version__ = "2022.03.14"



@define_app
def alignment_filter(align: c3_types.AlignedSeqsType, threshold: float=1.5) -> c3_types.AlignedSeqsType:
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


@define_app
def sequence_filter(seqs: c3_types.UnalignedSeqsType, threshold:float=50):
    """filters an alignment with a sequence with lots of repeats

    Parameters
    ----------
    seqs : Alignment
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
    if not seqs:
        return seqs

    # compute entropies for each value of k per seq
    ks = [3, 6, 9, 12]
    cevs = []
    for seq in seqs.seqs:
        entropies = [seq.counts(motif_length=k).entropy for k in ks]
        std = np.std(entropies)
        mean = np.mean(entropies)
        cevs.append(std / mean)

    cev = np.mean(cevs)
    if std != 0 and cev >= threshold:
        return NotCompleted(
            "FAIL",
            "sequence_filter",
            f"threshold < cev {cev} ",
            source=seqs.info.source,
        )

    seqs.info["sequence_cev"] = cev
    return seqs
