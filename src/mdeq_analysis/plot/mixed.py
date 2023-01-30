from pathlib import Path

from plotly.io import full_figure_for_development
from plotly.subplots import make_subplots

from .nabla import histogram_nabla_diff, fig_comparing_jsd_delta_nabla, nabla_vs_delta_nabla_traces
from .quantiles import get_one_plot, make_smiles_fig


def mixed_smiled_hist(ape=True, nbins=30):
    col = "bootstrap_pval"
    cats = ["CDS", "Intron"] if ape else ["Dmel", "Dsim"]
    smiles = [make_smiles_fig(cat, col) for cat in cats]
    # turn off the legend for the left-most subplot
    for e in smiles[0].data:
        e.showlegend = False

    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.5, 0.5],
        row_heights=[0.5, 0.5],
        vertical_spacing=0.15,
        subplot_titles=[f"(<b>{n}</b>) {c}" for n, c in zip("abc", cats + [""])],
        specs=[
            [{"type": "scatter"}, {"type": "scatter"}],
            [{"type": "histogram", "colspan": 2}, None],
        ],
    )
    for col, smile in enumerate(smiles, start=1):
        fig.add_traces(smile.data, cols=col, rows=1)
        fig.add_annotation(smile.layout.annotations[0], col=col, row=1)

    attr = dict(
        width=700,
        height=700,
        margin=dict(l=60, r=10, t=25, b=60),
    )
    fig.update_layout(legend_title_text="Data Type", **attr)
    fig.update_xaxes(
        title_font_size=18, tickfont=dict(size=14), title_standoff=5, dtick=0.25
    )
    fig.update_yaxes(
        title_font_size=18, tickfont=dict(size=14), title_standoff=5, dtick=0.25
    )

    hist = histogram_nabla_diff(ape=ape, nbins=nbins)
    fig.add_traces(hist.data, cols=1, rows=2)
    fig.update_layout(
        yaxis3=dict(title=hist.layout.yaxis.title),
        xaxis3=dict(title=hist.layout.xaxis.title),
    )

    fig.update_xaxes(
        title_font_size=18, tickfont=dict(size=14), title_standoff=5, dtick=0.25
    )
    fig.update_yaxes(
        title_font_size=18, tickfont=dict(size=14), title_standoff=5, dtick=0.25
    )
    fig.update_xaxes(showgrid=True, col=1, row=1, title="Theoretical Quantiles")
    fig.update_xaxes(showgrid=True, col=2, row=1, title="Theoretical Quantiles")
    fig.update_yaxes(
        showgrid=True, col=1, row=1, title=r"$p-\text{value}(\text{bootstrap})$"
    )
    fig.update_yaxes(
        showgrid=True,
        col=1,
        row=2,
        title="Probability",
        title_font_size=18,
        tickfont=dict(size=14),
        title_standoff=5,
        dtick=0.05,
    )
    fig.add_vline(x=0, row=2, line_width=3, line_dash="dash", line_color="black")
    # address plotly bug, suppress MathJax warning box
    full_figure_for_development(fig, warn=False)
    return fig


def make_mixed_properties(pval_paths, align_path, nabla_path, conv_paths):
    fig = make_subplots(
        rows=2,
        cols=3,
        column_widths=[0.33, 0.33, 0.33],
        row_heights=[0.5, 0.5],
        horizontal_spacing=0.07,
        vertical_spacing=0.05,
        subplot_titles=[f"(<b>{n}</b>)" for n in "abc"],
        specs=[
            [{"type": "scatter", "rowspan": 2}, {"type": "violin", "rowspan": 1}, {"type": "scatter", "rowspan": 2}],
            [None, {"type": "violin", "rowspan": 1}, None],
        ],
    )

    # left subplot
    left = get_one_plot(pval_paths, None, None)
    fig.add_traces(left.data, cols=1, rows=1)

    # middle piece, mid is a dict
    mid, stat_names = nabla_vs_delta_nabla_traces(conv_paths)
    fig.add_traces(mid[stat_names[0]], cols=2, rows=1)
    fig.add_traces(mid[stat_names[1]], cols=2, rows=2)

    # right plot
    right = fig_comparing_jsd_delta_nabla(align_path, nabla_path, None, None)
    fig.add_traces(right.data, cols=3, rows=1)

    # Setting layout
    fig.update_layout(
        legend_traceorder="reversed",
        legend_title_text="Alignment Length",
        legend_title_font_size=12,
        legend=dict(yanchor="top", y=0.98, xanchor="left", x=0.01),
    )
    shift = -35
    ## left
    fig.update_yaxes(range=[-0.05, 1.05], col=1, row=1, dtick=0.2)
    # y-axis
    fig.add_annotation(x=-0.05, y=0.5, yanchor="middle", xshift=shift, yref="y1", xref="x1", textangle=-90, showarrow=False,
                       text=r"$\large p-\text{value}(\chi^2)$")
    # x-axis
    fig.update_xaxes(range=[-0.05, 1.05], col=1, row=1, dtick=0.2)
    fig.add_annotation(x=0.5, y=-0.05, xanchor="center", yshift=shift, yref="y1", xref="x1", showarrow=False,
                       text=r"$\large \text{Theoretical Quantiles}$")

    # middle plots
    for row in (1, 2):
        fig.update_yaxes(range=[-0.3, 0.65], col=2, row=row)
        fig.update_xaxes(range=[-.5, 2.75], col=2, row=row, showticklabels=row==2)

    # y-axis titles
    # top
    fig.add_annotation(y=sum([-.3, .65]) / 2, x=-.5, yanchor="middle", xshift=shift, yref="y2", xref="x2", textangle=-90,
                       showarrow=False, text=r"$\large \hat{\nabla}$")
    # bottom
    fig.add_annotation(y=sum([-.3, .65]) / 2, x=-.5, yanchor="middle", xshift=shift, yref="y4", xref="x4", textangle=-90,
                       showarrow=False, text=r"$\large \hat\delta_{\nabla}$")
    # x-axis title bottom only
    fig.add_annotation(y=-0.3, x=sum([-.5, 2.75]) / 2, xanchor="center", yshift=shift, yref="y4", xref="x4",
                       showarrow=False, text=r"$\large \text{Alignment Length}$")

    # right axis titles
    fig.update_yaxes(range=[-0.1, 0.65], col=3, row=1)
    fig.update_xaxes(range=[-0.005, 0.06], col=3, row=1)
    fig.add_annotation(x=-0.005, y=0.3, yanchor="middle", xshift=shift, yref="y3", xref="x3", textangle=-90, showarrow=False,
                       text=r"$\large \hat\delta_{\nabla}$")
    fig.add_annotation(x=sum([-0.005, 0.06])/2, y=-.1, xanchor="center", yshift=shift, yref="y3", xref="x3", showarrow=False,
                       text=r"$\large \widehat{JSD}$")
    fig.update_xaxes(tickfont=dict(size=10))
    fig.update_yaxes(tickfont=dict(size=10))
    attr = dict(
        width=900,
        height=300,
        margin=dict(l=60, r=10, t=25, b=60),
    )
    fig.update_layout(**attr)
    return fig
