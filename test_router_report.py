import scanpy as sc

from backend.analysis_router import AnalysisRouter
from backend.router_report import RouterReport


DATA_FILE = "data/covid_pbmc/covid_pbmc.h5ad"

print("\nLoading dataset...")
adata = sc.read_h5ad(DATA_FILE)

router = AnalysisRouter()
profile = router.inspect(adata)
results = router.run(adata, profile=profile)

reporter = RouterReport()
report = reporter.build(results)

print("\n" + report[:8000])  # keep console output manageable

reporter.save(report, "reports/router_report.txt")
