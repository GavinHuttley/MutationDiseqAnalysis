from collections import defaultdict

from cogent3 import make_table
from cogent3.app import io
from cogent3.app import result as c3result
from mdeq.toe import ALT_TOE, NULL_TOE
from mdeq.utils import rich_display
from rich.progress import track


def get_pvalues(dstore) -> defaultdict:
    """returns chisq_pvals and bootstrap_pvals"""
    loader = io.load_db()
    pvals = defaultdict(list)
    is_hyp_result = None
    for m in track(dstore):
        obj = loader(m)
        is_hyp_result = isinstance(obj, c3result.hypothesis_result)
        if is_hyp_result:
            pval = obj.pvalue
        else:
            pval = obj.observed.get_hypothesis_result(NULL_TOE, ALT_TOE).pvalue

        if pval is None:
            # fitting issue
            continue

        pvals["chisq_pvals"].append(pval)
        if not is_hyp_result:
            pvals["bootstrap_pvals"].append(obj.pvalue)
        pvals["source"].append(obj.source)
    return pvals


def write_quantiles(path, limit=None, overwrite=False, verbose=0):
    """writes tsv file corresponding to path"""
    tsv_path = path.parent / f"{path.stem}.tsv"
    if tsv_path.exists() and not overwrite:
        return tsv_path

    if verbose:
        print(f"loading {limit} dstore records")

    limit = limit or 1000
    dstore = io.get_data_store(path, limit=limit)
    if len(dstore.incomplete) > 0:
        print(f"{path.stem} has {limit - len(dstore.incomplete)} incompleted results!")
        rich_display(dstore.summary_incomplete)
        if len(dstore) == 0:
            return False

    if verbose:
        print("getting p-values")

    pvals = get_pvalues(dstore)
    table = make_table(data=pvals)
    table.write(tsv_path)
    return path
