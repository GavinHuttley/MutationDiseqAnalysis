from collections import defaultdict
from pathlib import Path

from cogent3 import load_table
from cogent3.maths.stats.distribution import theoretical_quantiles
from cogent3.util.table import Table
from mdeq.utils import estimate_freq_null
from mdeq_analysis.plot import util
from plotly.graph_objects import Figure
from plotly.io import full_figure_for_development
from plotly.subplots import make_subplots


def load_quantiles(path, col="chisq_pval"):
    path = Path(path)
    tsv_path = path.parent / f"{path.stem}.tsv"

    table = load_table(tsv_path)
    table = table.get_columns(["name", col])

    for c, data in table.columns.items():
        if data.dtype == float:
            table.columns[c] = sorted(data)
    table.title = path.stem
    table.columns["theoretical"] = theoretical_quantiles(table.shape[0], dist="uniform")
    return table[:, [c for c in table.header if c != "name"]]


def get_trace(table: Table, col: str, name: str, alpha: float, _show_legend) -> dict:
    """get's scatter trace of quantiles"""
    color = util.get_colour_for_name(name, alpha)
    showlegend = name not in _show_legend
    _show_legend.add(name)
    return dict(
        x=sorted(table.columns["theoretical"]),
        y=sorted(table.columns[col]),
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
    show_legend = set()
    for path in paths:
        seed, bp = util.path_components(path)
        table = load_quantiles(path, col=stat)
        grouped_traces[seed].append(get_trace(table, stat, bp, alpha, show_legend))

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
        width=900,
        height=800,
        legend_traceorder="reversed",
        legend_title_text="<b>Alignment Length</b>",
        legend_font=dict(size=14),
        margin=dict(l=60, r=10, t=25, b=60),
    )
    fig.update_xaxes(tickfont=dict(size=14), dtick=0.25)
    fig.update_yaxes(tickfont=dict(size=14), dtick=0.25)
    # address plotly bug, suppress MathJax warning box
    full_figure_for_development(fig, warn=False)
    return fig


def get_one_plot(paths, width, height):
    show_legend = set()
    d = {p.stem.split("-")[-2]: load_quantiles(p) for p in paths}
    traces = [get_trace(t, "chisq_pval", k, 0.8, show_legend) for k, t in d.items()]
    for t in traces:
        t["showlegend"] = True

    fig = Figure(data=traces)
    fig.update_layout(
        width=width,
        height=height,
        legend_traceorder="reversed",
        xaxis_title=r"$\text{Theoretical Quantiles}$",
        yaxis_title=r"$p-\text{value}(\chi^2)$",
        legend_title_text="Alignment Length",
        legend_title_font_size=16,
        margin=dict(l=60, r=10, t=25, b=25),
        legend=dict(yanchor="top", y=0.97, xanchor="left", x=0.03),
    )
    # address plotly bug, suppress MathJax warning box
    full_figure_for_development(fig, warn=False)
    fig.update_xaxes(title_font_size=18, title_standoff=5)
    fig.update_yaxes(title_font_size=18, title_standoff=5)
    return fig


def load_dros_quantiles(spec: str, col: str) -> tuple[Table, Table, Table]:
    """loads probability quantiles from file

    Parameters
    ----------
    spec : str
        Drosophila species abbreviation, Dmel or Dsim
    col : str
        table column containing probabilities of interest

    Returns
    -------
    tuple[Table, Table, Table]
        observed, negative control, positive control
    """
    obs = load_quantiles(
        f"../results/drosophila/fg-GSN-toe/toe-dmel_dsim_dyak-{spec}.tsv", col=col
    )
    pos = load_quantiles(
        f"../results/drosophila/toe-controls-results/toe-pos_control-toe-dmel_dsim_dyak-{spec}.tsv",
        col=col,
    )
    neg = load_quantiles(
        f"../results/drosophila/toe-controls-results/toe-neg_control-toe-dmel_dsim_dyak-{spec}.tsv",
        col=col,
    )
    return obs, neg, pos


def load_ape_quantiles(klass: str, col: str) -> tuple[Table, Table, Table]:
    """loads probability quantiles from file

    Parameters
    ----------
    klass : str
        sequence class, cds or intron
    col : str
        table column containing probabilities of interest

    Returns
    -------
    tuple[Table, Table, Table]
        observed, negative control, positive control
    """
    obs = load_quantiles(
        f"../results/ape/toe/fg-GSN-toe/filtered-ape_aligned_{klass}.tsv", col=col
    )
    pos = load_quantiles(
        f"../results/ape/toe-controls-results/toe-pos_control-filtered-ape_aligned_{klass}.tsv",
        col=col,
    )
    neg = load_quantiles(
        f"../results/ape/toe-controls-results/toe-neg_control-filtered-ape_aligned_{klass}.tsv",
        col=col,
    )
    return obs, neg, pos


def make_smiles_fig(cat: str, col: str) -> Figure:
    """makes scatter plot showing the negative control, observed and positive control quantiles for th data indicated by cat

    Parameters
    ----------
    cat : str
        one of Dmel, Dsim, cds, intron
    col : str
        table column name, either chisq_pval or bootstrap_pval

    Returns
    -------
    Figure
    """
    loader = load_dros_quantiles if cat in "DsimDmel" else load_ape_quantiles
    obs, neg, pos = loader(cat, col)

    show_legend = set()
    obs_trace = get_trace(obs, col, "Observed", 0.8, show_legend)
    # the freq of neutral mutation diseq, f_NMD
    f_MD = 1 - estimate_freq_null(obs.columns["bootstrap_pval"])
    pos_trace = get_trace(pos, col, "+ve", 0.8, show_legend)
    neg_trace = get_trace(neg, col, "-ve", 0.8, show_legend)
    fig = Figure(
        {
            "data": [neg_trace, obs_trace, pos_trace],
            "layout": {"title": cat, "width": 600, "height": 600, "showlegend": True},
        }
    )
    latex = r"\hat f_{\text{NMD}}\approx"
    fig.add_annotation(
        text=f"${latex}{f_MD:.2f}$",
        x=0.1,
        y=0.8,
        xref="x domain",
        yref="y domain",
        showarrow=False,
        font=dict(size=20),
    )
    return fig


def subplot_smiles(ape=True):
    col = "bootstrap_pval"
    cats = ["CDS", "Intron"] if ape else ["Dmel", "Dsim"]
    smiles = [make_smiles_fig(cat, col) for cat in cats]
    # turn off the legend for the left-most subplot
    for e in smiles[0].data:
        e.showlegend = False

    fig = make_subplots(
        rows=1,
        cols=2,
        shared_xaxes=True,
        shared_yaxes=True,
        vertical_spacing=0.05,
        horizontal_spacing=0.05,
        x_title="Theoretical Quantiles",
        y_title=r"$p-\text{value (bootstrap})$",
        subplot_titles=cats,
    )
    for col, smile in enumerate(smiles, start=1):
        fig.add_traces(smile.data, cols=col, rows=1)

    attr = dict(
        width=700,
        height=400,
        legend_title_text="Data Type",
        margin=dict(l=60, r=10, t=25, b=60),
    )
    fig.update_layout(**attr)
    fig.update_xaxes(
        title_font_size=18, tickfont=dict(size=14), title_standoff=5, dtick=0.25
    )
    fig.update_yaxes(
        title_font_size=18, tickfont=dict(size=14), title_standoff=5, dtick=0.25
    )
    # address plotly bug, suppress MathJax warning box
    full_figure_for_development(fig, warn=False)
    return fig
