import os
from pathlib import Path

import pandas as pd
import scanpy as sc

from backend.gene_mapper import GeneMapper
from backend.condition_de import ConditionDE
from backend.condition_pathways import ConditionPathways
from backend.condition_interpreter import ConditionInterpreter
from backend.volcano_plot import VolcanoPlot
from backend.de_heatmap import DEHeatmap


DATA_FILE = "data/covid_pbmc/covid_pbmc.h5ad"
RANKING_FILE = "reports/celltype_covid_response_ranking.csv"
OUT_DIR = Path("outputs")
REPORT_DIR = Path("reports")
REPORT_FILE = REPORT_DIR / "covid_top3_driver_report.txt"

OUT_DIR.mkdir(exist_ok=True, parents=True)
REPORT_DIR.mkdir(exist_ok=True, parents=True)


def safe_name(text: str) -> str:
    return (
        str(text)
        .replace("/", "_")
        .replace(" ", "_")
        .replace(".", "_")
        .replace(":", "_")
    )


def load_dataset() -> sc.AnnData:
    print("\nLoading real COVID dataset...")
    adata = sc.read_h5ad(DATA_FILE)

    print("\nDataset loaded.")
    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")

    mapper = GeneMapper()
    adata = mapper.fix_gene_names(adata)

    if "celltype.final" not in adata.obs.columns:
        raise ValueError("Expected 'celltype.final' in adata.obs.")

    if "disease" not in adata.obs.columns:
        raise ValueError("Expected 'disease' in adata.obs.")

    adata.obs["cell_type"] = adata.obs["celltype.final"].astype(str)
    adata.obs["condition"] = adata.obs["disease"].astype(str)

    adata = adata[adata.obs["condition"].isin(["normal", "COVID-19"])].copy()
    adata.obs["condition"] = adata.obs["condition"].replace(
        {"normal": "Healthy", "COVID-19": "COVID"}
    )

    print("\nCondition counts:")
    print(adata.obs["condition"].value_counts())

    print("\nCell type counts:")
    print(adata.obs["cell_type"].value_counts().head(15))

    return adata


def build_or_load_ranking(adata: sc.AnnData) -> pd.DataFrame:
    if os.path.exists(RANKING_FILE):
        print(f"\nLoading existing ranking file: {RANKING_FILE}")
        ranking = pd.read_csv(RANKING_FILE)
        ranking = ranking.sort_values("de_genes", ascending=False).reset_index(drop=True)
        return ranking

    print("\nRanking file not found. Building ranking from scratch...")

    de_engine = ConditionDE()
    pathway_engine = ConditionPathways()
    interpreter = ConditionInterpreter()

    results = []

    celltypes = (
        adata.obs["cell_type"]
        .value_counts()
    )
    celltypes = celltypes[celltypes >= 200].index.tolist()

    for celltype in celltypes:
        print(f"\nAnalyzing {celltype}")

        subset = adata[adata.obs["cell_type"] == celltype].copy()
        counts = subset.obs["condition"].value_counts()

        if counts.get("COVID", 0) < 50 or counts.get("Healthy", 0) < 50:
            continue

        subset = de_engine.run(
            subset,
            condition_key="condition",
            group1="COVID",
            group2="Healthy"
        )

        de = de_engine.top_genes(
            subset,
            group="COVID",
            n_genes=200
        )

        sig = de[
            (de["pvals_adj"] < 0.05) &
            (de["logfoldchanges"].abs() > 0.5)
        ].copy()

        pathways = pathway_engine.analyze(sig, top_n=50)
        interpretation = interpreter.interpret(pathways)

        results.append({
            "cell_type": celltype,
            "de_genes": 0 if sig.empty else len(sig),
            "covid_cells": int(counts["COVID"]),
            "healthy_cells": int(counts["Healthy"]),
            "interpretation": "; ".join(interpretation),
        })

    ranking = pd.DataFrame(results)
    ranking = ranking.sort_values("de_genes", ascending=False).reset_index(drop=True)
    ranking.to_csv(RANKING_FILE, index=False)

    print(f"\nSaved ranking file: {RANKING_FILE}")
    return ranking


def run_top3_analysis(adata: sc.AnnData, ranking: pd.DataFrame) -> list[dict]:
    de_engine = ConditionDE()
    pathway_engine = ConditionPathways()
    interpreter = ConditionInterpreter()
    volcano = VolcanoPlot()
    heatmap = DEHeatmap()

    top3 = ranking.head(3).copy()

    summaries = []

    for _, row in top3.iterrows():
        celltype = row["cell_type"]
        print(f"\nRunning deep analysis for top driver: {celltype}")

        subset = adata[adata.obs["cell_type"] == celltype].copy()
        counts = subset.obs["condition"].value_counts()

        if counts.get("COVID", 0) < 50 or counts.get("Healthy", 0) < 50:
            print(f"Skipping {celltype}: not enough cells in one condition.")
            continue

        subset = de_engine.run(
            subset,
            condition_key="condition",
            group1="COVID",
            group2="Healthy"
        )

        de = de_engine.top_genes(
            subset,
            group="COVID",
            n_genes=200
        )

        sig = de[
            (de["pvals_adj"] < 0.05) &
            (de["logfoldchanges"].abs() > 0.5)
        ].copy()

        if sig.empty:
            sig_for_pathways = de.head(50).copy()
        else:
            sig_for_pathways = sig.copy()

        pathways = pathway_engine.analyze(
            sig_for_pathways,
            top_n=50
        )

        if pathways is not None and not pathways.empty:
            pathway_terms = pathways["Term"].head(10).tolist()
            interpretation = interpreter.interpret(pathways)
        else:
            pathway_terms = []
            interpretation = ["No strong condition-specific pathway signature identified."]

        safe = safe_name(celltype)

        if de is not None and not de.empty:
            volcano.create(
                de,
                output_file=str(OUT_DIR / f"volcano_{safe}_COVID_vs_Healthy.png"),
                title=f"Volcano plot: {celltype} (COVID vs Healthy)"
            )

        top_heatmap_genes = (
            sig_for_pathways["names"].head(15).tolist()
            if sig_for_pathways is not None and not sig_for_pathways.empty
            else de["names"].head(15).tolist()
        )

        if top_heatmap_genes:
            heatmap.create(
                subset,
                top_heatmap_genes,
                groupby="condition",
                output_file=str(OUT_DIR / f"heatmap_{safe}_COVID_vs_Healthy.png")
            )

        summaries.append({
            "cell_type": celltype,
            "covid_cells": int(counts["COVID"]),
            "healthy_cells": int(counts["Healthy"]),
            "de_genes": 0 if sig.empty else len(sig),
            "top_genes": sig_for_pathways["names"].head(10).tolist() if sig_for_pathways is not None and not sig_for_pathways.empty else [],
            "top_pathways": pathway_terms,
            "interpretation": interpretation,
        })

    return summaries


def write_report(ranking: pd.DataFrame, summaries: list[dict]) -> None:
    lines = []
    lines.append("========== OMNISEQAI COVID DRIVER REPORT ==========\n")
    lines.append(f"Generated from: {DATA_FILE}")
    lines.append(f"Cells analyzed: {sum(r['covid_cells'] + r['healthy_cells'] for r in summaries):,}")
    lines.append("")

    lines.append("===== TOP COVID-RESPONSIVE CELL TYPES =====")
    lines.append(ranking.head(10).to_string(index=False))

    for item in summaries:
        lines.append("\n" + "=" * 80)
        lines.append(f"Cell type: {item['cell_type']}")
        lines.append(f"COVID cells: {item['covid_cells']}")
        lines.append(f"Healthy cells: {item['healthy_cells']}")
        lines.append(f"Significant DE genes: {item['de_genes']}")

        lines.append("\nTop DE genes:")
        if item["top_genes"]:
            lines.append(", ".join(item["top_genes"]))
        else:
            lines.append("None")

        lines.append("\nTop pathways:")
        if item["top_pathways"]:
            for p in item["top_pathways"]:
                lines.append(f"- {p}")
        else:
            lines.append("- None")

        lines.append("\nInterpretation:")
        for entry in item["interpretation"]:
            lines.append(f"- {entry}")

    with open(REPORT_FILE, "w") as f:
        f.write("\n".join(lines))

    print(f"\nReport saved: {REPORT_FILE}")


def main():
    adata = load_dataset()
    ranking = build_or_load_ranking(adata)

    if ranking.empty:
        print("\nNo ranking results found. Exiting.")
        return

    summaries = run_top3_analysis(adata, ranking)
    write_report(ranking, summaries)

    print("\n===== TOP 3 DRIVER SUMMARY =====")
    for item in summaries:
        print("\n" + "-" * 70)
        print(f"Cell type: {item['cell_type']}")
        print(f"Significant DE genes: {item['de_genes']}")
        print("Interpretation:")
        for entry in item["interpretation"]:
            print(f"- {entry}")

    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
