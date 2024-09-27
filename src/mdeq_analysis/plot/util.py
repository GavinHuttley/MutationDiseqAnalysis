import re
import time
from plotly.io import write_image
from pathlib import Path

from cogent3.maths.measure import jsd

_seed = re.compile("(hi|lo)_(hi|lo)")
_aln_length = re.compile(r"\d+bp")


def path_components(path: Path) -> tuple[str]:
    """returns the JSD/Entrop label and alignment length from file path"""
    bp = _aln_length.findall(path.stem)[0]
    seed = "_".join(_seed.findall(path.stem)[0])
    return seed, bp


def keyed_traces(traces) -> dict:
    """keyed by trace name"""
    return {x["name"]: x for x in traces}


def annotated_subplot(fig, seed, col, row, x, y):
    """puts text for hi/low entropy on subplot component of fig"""
    text = {
        "lo_lo": "Low JSD, Low Entropy",
        "lo_hi": "Low JSD, High Entropy",
        "hi_lo": "High JSD, Low Entropy",
        "hi_hi": "High JSD, High Entropy",
    }[seed]
    fig.add_annotation(
        text=text,
        row=row,
        col=col,
        xref="x domain",
        yref="y domain",
        x=x,
        y=y,
        showarrow=False,
        font=dict(size=14),
    )


def get_colour_for_name(name: str, alpha: float = 0.4) -> str:
    """returns color for length"""
    return {
        "300bp": f"rgba(0, 128, 255, {alpha})",
        "3000bp": f"rgba(253, 214, 3,{alpha})",
        "30000bp": f"rgba(209,17,65,{alpha})",
        "+ve": f"rgba(214,45,32,{alpha})",
        "-ve": f"rgba(0,87,231,{alpha})",
        "Observed": f"rgba(0,0,0,{alpha})",
    }[name]


def calc_jsd(aln):
    counts = aln.counts_per_seq()
    freqs = counts.to_freq_array()
    return jsd(*freqs.array)


class pdf_writer:
    """class that handles super annooying mathjax warning box in plotly pdf's"""

    def __init__(self) -> None:
        self._done_once = False

    def __call__(self, fig, path):
        # the sleep, plus successive write, is ESSENTIAL to avoid the super annoying
        # "[MathJax]/extensions/MathMenu.js" text box error
        # but we only need to do this once
        if not self._done_once:
            write_image(fig, path)
            time.sleep(2)
            self._done_once = True
        
        write_image(fig, path)
