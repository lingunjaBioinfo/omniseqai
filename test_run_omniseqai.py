import subprocess
import sys


DATA_FILE = "data/covid_pbmc/covid_pbmc.h5ad"


subprocess.run(
    [
        sys.executable,
        "run_omniseqai.py",
        "--input",
        DATA_FILE,
        "--output",
        "reports/router_report.txt",
        "--preview-lines",
        "120"
    ],
    check=True
)
