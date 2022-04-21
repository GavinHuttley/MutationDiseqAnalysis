from collections import defaultdict
from pathlib import Path

import plotly.express as px

from cogent3 import make_table
from cogent3.app import io
from mdeq.sqlite_data_store import sql_loader
from plotly.subplots import make_subplots

from mdeq_analysis.plot import util


def convert_to_table(path):
    """converts the delta_nabla instances to a table"""
    from mdeq import convergence  # required to register the deserialiser

    loader = sql_loader()
    dstore = io.get_data_store(path)
    results = defaultdict(list)
    for m in dstore:
        result = loader(m)
        results["nabla"].append(result.obs_nabla)
        results["delta_nabla"].append(result.delta_nabla)
    return make_table(data=results, title=path.stem)


def stat_to_trace(table, col, name):
    """returns violin plot trace"""
    color = util.get_colour_for_name(name)
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
        grouped_traces[seed].append(stat_to_trace(table, stat, bp))

    y_title = r"$\hat\nabla$" if stat == "nabla" else r"$\hat\delta_{\nabla}$"
    title = (
        r"$\nabla \text{ standard scaling}$"
        if stat == "nabla"
        else r"$\delta_{nabla} \text{ standard scaling}$"
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
    return fig
