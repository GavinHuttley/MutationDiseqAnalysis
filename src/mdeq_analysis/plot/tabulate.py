import re

from cogent3 import load_table, make_table, open_data_store
from mdeq.utils import load_from_sqldb

_num = re.compile(r"\d+")


def get_pvalues(path):
    table = load_table(path)
    table.columns["rank"] = [_num.findall(v)[0] for v in table.columns["name"]]
    return table[:, ["rank", "bootstrap_pval"]]


def get_delta_nabla(path):
    dstore = open_data_store(path)
    loader = load_from_sqldb()
    rows = []
    for m in dstore:
        r = loader(m)
        rows.append((r.source, r.delta_nabla, r.std_null))
    table = make_table(header=["name", "delta_nabla", "std"], data=rows)
    table.columns["rank"] = [_num.findall(v)[0] for v in table.columns["name"]]
    return table[:, ["rank", "delta_nabla", "std"]]


def get_alignment_lengths(path):
    dstore = open_data_store(path)
    loader = load_from_sqldb()
    rows = []
    for m in dstore:
        aln = loader(m)
        mmu = aln.get_lengths()["mmu"]
        rows.append([m.unique_id, mmu])

    table = make_table(header=["name", "length"], data=rows)
    table.columns["rank"] = [_num.findall(v)[0] for v in table.columns["name"]]
    return table[:, ["rank", "length"]]


def merged(pvals, delta_nabla, align_lengths):
    pvals.index_name = delta_nabla.index_name = align_lengths.index_name = "rank"

    m = pvals.inner_join(delta_nabla)
    m = m.inner_join(align_lengths, digits=2)
    old_cols = [c for c in m.columns if c.startswith("right")]
    new_cols = [c.replace("right_", "") for c in old_cols]
    m = m.with_new_header(old_cols, new_cols)
    return m


def _format_element(x):
    return f"{x:.2f}" if type(x) == float else x


def make_latex_table(table):
    """changes to row"""
    table.columns["length"] = [f"{e:,}" for e in table.columns["length"]]
    m = table.transposed("", select_as_header="rank", digits=2)
    for c, col in m.columns.items():
        m.columns[c] = [_format_element(e) for e in col]

    data = m.columns.to_dict()
    data[r""] = [
        r"$\hat{p}$-value",
        r"$\hat\nabla_c$",
        r"$\hat\sigma_\nabla$",
        "length",
    ]
    table = make_table(data=data)
    new_order = [table.header[0]] + sorted(table.header[1:])
    table = table[:, new_order]
    return table.with_new_header(
        table.header, (r"Statistic \textbackslash ~ Intron rank",) + table.header[1:]
    )


def fxy_table(data_dir, result_dir):
    align_lengths = get_alignment_lengths(
        data_dir / "introns-aligned-filtered.sqlitedb"
    )
    pvalues = get_pvalues(result_dir / "toe/toe-fxy-intron-mmu.tsv")
    dnabla = load_table(
        result_dir / "convergence/convergence-fxy-intron-mmu.tsv"
    )
    dnabla.columns["rank"] = [_num.findall(v)[0] for v in dnabla.columns["source"]]
    dnabla = dnabla[:, ["rank", "delta_nabla", "std"]].sorted(columns="rank")
    return make_latex_table(merged(pvalues, dnabla, align_lengths))
