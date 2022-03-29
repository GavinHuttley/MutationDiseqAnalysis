import re

from pathlib import Path


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


def get_colour_for_length(length: str, alpha: float = 0.4) -> str:
    """returns color for length"""
    return {
        "300bp": f"rgba(0,174,219,{alpha})",
        "3000bp": f"rgba(255,196,37,{alpha})",
        "30000bp": f"rgba(209,17,65,{alpha})",
    }[length]
