from collections import defaultdict
from pathlib import Path

from cogent3 import load_table
from cogent3.maths.stats.distribution import theoretical_quantiles
from plotly.subplots import make_subplots

from . import util


def load_quantiles(path, col="chisq_pvals"):
    tsv_path = path.parent / f"{path.stem}.tsv"

    table = load_table(tsv_path)
    for col, data in table.columns.items():
        if data.dtype == float:
            table.columns[col] = sorted(data)
    table.title = path.stem
    table.columns["theoretical"] = theoretical_quantiles(table.shape[0], dist="uniform")
    return table[:, [c for c in table.header if c != "source"]]


_show_legend = set()


def get_trace(table, col, name, alpha):
    """get's scatter trace of quantiles"""
    color = util.get_colour_for_length(name, alpha)
    showlegend = name not in _show_legend
    _show_legend.add(name)
    return dict(
        x=table.columns["theoretical"],
        y=table.columns[col],
        name=name,
        type="scatter",
        mode="lines",
        line=dict(
            width=2,
            color=color,
            shape="spline",
            smoothing=1.3,
        ),
        marker=dict(color=color, size=4),
        legendgroup=name,
        showlegend=showlegend,
    )


def get_quantile_fig(paths: Path, stat: str, alpha: float = 0.4):
    """returns scatter subplot of quantiles"""
    grouped_traces = defaultdict(list)
    for path in paths:
        seed, bp = util.path_components(path)
        table = load_quantiles(path)
        grouped_traces[seed].append(get_trace(table, stat, bp, alpha))

    fig = make_subplots(
        rows=2,
        cols=2,
        shared_xaxes=True,
        shared_yaxes=True,
        x_title="Theoretical Quantiles",
        y_title=r"$p-\text{value }(\chi^2)$",
        vertical_spacing=0.05,
        horizontal_spacing=0.05,
    )

    # <JSD>_<Entropy> -> col, row (which is array coordinates)
    col_row = dict(lo_hi=(1, 1), hi_hi=(2, 1), lo_lo=(1, 2), hi_lo=(2, 2))

    for seed, (col, row) in col_row.items():
        traces = util.keyed_traces(grouped_traces[seed])
        # the following order is because of plotly behaviour
        # it produces the legend in the correct order
        for size in ("300bp", "30000bp", "3000bp"):
            trace = traces.get(size, None)
            if trace is None:
                continue
            fig.add_trace(trace, col=col, row=row)
        util.annotated_subplot(fig, seed, col, row, 0.1, 0.9)

    fig.update_layout(
        width=800,
        height=800,
        legend_traceorder="reversed",
        legend_title_text="<b>Alignment Length</b>",
        legend_font=dict(size=14),
    )
    fig.update_xaxes(tickfont=dict(size=14), dtick=0.25)
    fig.update_yaxes(tickfont=dict(size=14), dtick=0.25)
    return fig
