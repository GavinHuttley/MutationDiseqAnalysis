import pathlib

ROOT_DIR = pathlib.Path(__file__).parent.parent.parent

OUTPUT_ROOT = ROOT_DIR / "MutationDiseqMS"

FIG_DIR = OUTPUT_ROOT / "figs"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SUPP_FIG_DIR = OUTPUT_ROOT / "figs_supp"
SUPP_FIG_DIR.mkdir(parents=True, exist_ok=True)

TABLE_DIR = OUTPUT_ROOT / "tables"
TABLE_DIR.mkdir(parents=True, exist_ok=True)

SUPP_TABLE_DIR = OUTPUT_ROOT / "tables_supp"
SUPP_TABLE_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = ROOT_DIR / "MutationDiseqAnalysis/data"
RESULT_DIR = ROOT_DIR / "MutationDiseqAnalysis/results"
