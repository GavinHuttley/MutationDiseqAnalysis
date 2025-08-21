# download the data and result files
import pathlib
import shutil
import urllib.request
import zipfile

root_dir = pathlib.Path(__file__).parent
while not (root_dir / "MutationDiseqAnalysis").exists():
    root_dir = root_dir.parent

root_dir = root_dir / "MutationDiseqAnalysis"

DATA_URL = "https://zenodo.org/records/16916309/files/mdeq_data.zip?download=1"
DATA_NAME = "mdeq_data.zip"
RESULT_URL = "https://zenodo.org/records/16916309/files/mdeq_results.zip?download=1"
RESULT_NAME = "mdeq_results.zip"


def get_install_remote(url: str, dest_zip: str, dest: str) -> str:
    zip_dest = root_dir / dest_zip
    unzipped_dest = root_dir / zip_dest.stem
    expected = root_dir / dest
    if expected.exists():
        return dest

    if zip_dest.exists() and expected.exists():
        # we will inflate zip archive each time
        shutil.rmtree(expected)
    elif not zip_dest.exists():
        urllib.request.urlretrieve(url, filename=zip_dest)  # noqa: S310

    with zipfile.ZipFile(zip_dest, "r") as zip_ref:
        zip_ref.extractall(root_dir)

    unzipped_dest.rename(expected)
    return dest


if __name__ == "__main__":
    # get the data
    get_install_remote(DATA_URL, DATA_NAME, "data")
    # get the results -- this is ~20GB!
    get_install_remote(RESULT_URL, RESULT_NAME, "results")
