import os
import pathlib

if os.environ["LOGNAME"] == "gavin":
    ROOT_DIR = pathlib.Path("~/repos/MutationDiseq").expanduser()
else:
    ROOT_DIR = pathlib.Path("~/MutDiseq").expanduser()

OUTPUT_ROOT = ROOT_DIR / "MutationDiseqMS"

FIG_DIR = OUTPUT_ROOT / "figs"
FIG_DIR.mkdir(parents=True, exist_ok=True)

TABLE_DIR = OUTPUT_ROOT / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

SUPP_TABLE_DIR = OUTPUT_ROOT / "tables_supp"
SUPP_TABLE_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = ROOT_DIR / "MutationDiseqAnalysis/data"
RESULT_DIR = ROOT_DIR / "MutationDiseqAnalysis/results"
