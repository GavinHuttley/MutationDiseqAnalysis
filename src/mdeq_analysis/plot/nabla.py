from collections import defaultdict
from pathlib import Path

import plotly.express as px
from cogent3 import make_table, open_data_store
from mdeq.utils import load_from_sqldb
from mdeq_analysis.plot import util
from plotly.io import full_figure_for_development
from plotly.subplots import make_subplots


def convert_to_table(path):
    """converts the delta_nabla instances to a table"""
    # need import to register deserialiser
    from mdeq import convergence  # noqa: F401

    loader = load_from_sqldb()
    dstore = open_data_store(path)
    results = defaultdict(list)
    for m in dstore.completed:
        result = loader(m)
        results["nabla"].append(result.obs_nabla)
        results["nabla_c"].append(result.nabla_c)
    return make_table(data=results, title=path.stem)


def stat_to_trace(table, col, name, alpha=0.4):
    """returns violin plot trace"""
    color = util.get_colour_for_name(name, alpha=alpha)
    return dict(
        type="violin",
        y=table.columns[col],
        showlegend=False,
        name=name,
        fillcolor=color,
        line_color=color,
    )


def get_fig_for_stat(paths, stat):
    """returns violin subplots of convergence stat"""
    grouped_traces = defaultdict(list)
    for path in paths:
        table = convert_to_table(path)
        seed, bp = util.path_components(path)
        grouped_traces[seed].append(stat_to_trace(table, stat, bp, alpha=0.8))

    y_title = r"$\hat\nabla$" if stat == "nabla" else r"$\hat\nabla_c$"
    title = (
        r"$\nabla \text{ standard scaling}$"
        if stat == "nabla"
        else r"$\hat\nabla_c \text{ standard scaling}$"
    )
    fig = make_subplots(
        rows=2,
        cols=2,
        shared_xaxes=True,
        shared_yaxes=True,
        x_title="Alignment Length",
        y_title=y_title,
        vertical_spacing=0.05,
        horizontal_spacing=0.05,
    )
    col_row = dict(lo_hi=(1, 1), hi_hi=(2, 1), lo_lo=(1, 2), hi_lo=(2, 2))
    for seed, (col, row) in col_row.items():
        traces = util.keyed_traces(grouped_traces[seed])
        for size in ("300bp", "3000bp", "30000bp"):
            trace = traces.get(size, None)
            if trace is None:
                continue
            fig.add_trace(trace, col=col, row=row)
        util.annotated_subplot(fig, seed, col, row, 0.9, 0.9)

    fig.update_layout(
        title=title,
        width=800,
        height=800,
    )
    fig.update_xaxes(tickfont=dict(size=14), dtick=0.25)
    fig.update_yaxes(tickfont=dict(size=14), dtick=0.25)
    fig.layout.annotations[1].xshift = -50
    for annot in fig.layout.annotations:
        if "Entropy" in annot["text"]:
            continue
        annot["font"]["size"] = 18
    # address plotly bug, suppress MathJax warning box
    full_figure_for_development(fig, warn=False)
    return fig


def fig_nabla_vs_delta_nabla(paths, width, height):
    """groups results by statistic"""
    grouped_traces, stats = nabla_vs_delta_nabla_traces(paths)

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        shared_yaxes=False,
        vertical_spacing=0.05,
    )
    fig.add_traces(grouped_traces[stats[0]], cols=1, rows=1)
    fig.add_traces(grouped_traces[stats[1]], cols=1, rows=2)
    fig.update_xaxes(
        tickfont=dict(size=14),
        dtick=0.25,
        title_text=r"$\text{Alignment Length}$",
        title_font_size=18,
        title_standoff=5,
        row=2,
        col=1,
    )
    fig.update_yaxes(
        tickfont=dict(size=14),
        dtick=0.25,
        title_font_size=18,
        title_standoff=5,
    )
    fig.update_yaxes(title_text=r"$\hat{\nabla}$", col=1, row=1, title_font_size=18)
    fig.update_yaxes(
        title_text=r"$\hat\nabla_c$",
        col=1,
        row=2,
    )
    fig.update_layout(width=width, height=height, margin=dict(l=25, r=20, t=25, b=45))

    # fig.layout.annotations[-1].yshift = -20
    full_figure_for_development(fig, warn=False)
    return fig


def nabla_vs_delta_nabla_traces(paths):
    grouped_traces = defaultdict(list)
    stats = ("nabla", "nabla_c")
    for size in ("300bp", "3000bp", "30000bp"):
        for path in paths:
            if size in path.name:
                break
        else:
            raise ValueError(f"no path found for {size} in {paths}")
        table = convert_to_table(path)
        _, bp = util.path_components(path)
        for stat in stats:
            grouped_traces[stat].append(stat_to_trace(table, stat, bp, alpha=0.8))
    return grouped_traces, stats


def fig_comparing_jsd_delta_nabla(
    align_path, nabla_path, width, height, log_scale_x=False
):
    from .util import calc_jsd

    aligns = open_data_store(align_path)
    nablas = open_data_store(nabla_path)
    loader = load_from_sqldb()
    x, y = [], []
    for m in aligns.completed:
        aln = loader(m)
        n_m = [r for r in nablas.completed if r.unique_id == m.unique_id][0]
        nabla = loader(n_m)
        x.append(calc_jsd(aln))
        y.append(nabla.nabla_c)

    fig = px.scatter(x=x, y=y, opacity=0.7)
    fig.update_traces(marker={"size": 6})
    fig.update_layout(
        width=width,
        height=height,
        margin=dict(l=20, r=20, t=25, b=25),
    )
    xaxis_kwargs = dict(
        title_text=r"$\widehat{JSD}$",
        title_font_size=18,
        tickfont=dict(size=14),
        title_standoff=5,
    )
    if log_scale_x:
        xaxis_kwargs["type"] = "log"
    fig.update_xaxes(**xaxis_kwargs)
    fig.update_yaxes(
        title_text=r"$\hat\nabla_c$",
        title_font_size=18,
        tickfont=dict(size=14),
        title_standoff=5,
    )
    full_figure_for_development(fig, warn=False)
    return fig


def compare_nabla(ape=True):
    if ape:
        path_template = (
            "../results/ape/convergence/convergence-filtered-ape_aligned_{}.sqlitedb"
        )
        cats = ["CDS", "Intron"]
    else:
        path_template = "../results/drosophila/convergence/convergence-toe-dmel_dsim_dyak-{}.sqlitedb"
        cats = ["Dmel", "Dsim"]

    loader = load_from_sqldb()
    tables = []
    for cat in cats:
        rows = []
        dstore = open_data_store(path_template.format(cat))
        for m in dstore.completed:
            r = loader(m)
            name = Path(r.source).name.split(".")[0]
            rows.append((name, r.nabla_c))
        tables.append(make_table(["name", cat], data=rows, index_name="name"))

    table = tables[0].inner_join(tables[1])
    table = table.with_new_header(f"right_{cats[1]}", cats[1])
    table.columns["diff"] = table.columns[cats[0]] - table.columns[cats[1]]
    num_gt = sum(table.columns["diff"] > 0)
    num_gt = sum(table.columns["diff"] > 0)
    print(f"Number {cats[0]} > {cats[1]} = {num_gt} out of total {table.shape[0]}")
    return table


def histogram_nabla_diff(ape=True, nbins=30):
    table = compare_nabla(ape=ape)
    stat = r"\hat\nabla_c"
    elements = [r"\text{%s}" % c for c in table.columns if c not in "namediff"]
    for index in (0, -1):
        elements.insert(index, stat)

    axis_title = "${}({})-{}({})$".format(*elements)
    fig = px.histogram(x=table.columns["diff"], histnorm="probability", nbins=nbins)
    common_kwargs = dict(
        title_font_size=18, tickfont=dict(size=14), title_standoff=5, dtick=0.25
    )
    fig.update_xaxes(title=axis_title, **common_kwargs)
    fig.update_yaxes(title="Probability", **common_kwargs)
    attr = dict(
        width=700,
        height=400,
        margin=dict(l=60, r=10, t=25, b=60),
    )
    fig.update_layout(**attr)
    # fig.add_trace(px.scatter(x=[0,0], y=[1, 1]))
    return fig
