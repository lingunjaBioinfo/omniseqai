from backend.final_report import FinalReport
from backend.pipeline import SingleCellPipeline
from backend.annotation import CellAnnotator
from backend.markers import MarkerAnalyzer
from backend.pathway_analysis import PathwayAnalyzer
from backend.disease_interpreter import DiseaseInterpreter
from backend.conclusion_engine import ConclusionEngine
from backend.cell_communication import CellCommunication


# -------------------
# PIPELINE
# -------------------

pipeline = SingleCellPipeline()

adata = pipeline.load_pbmc3k()

pipeline.calculate_qc()

pipeline.filter_data()

pipeline.normalize_data()

pipeline.identify_hvg()

pipeline.run_pca()

pipeline.compute_neighbors()

pipeline.run_umap()

pipeline.cluster_cells()

# Use the fully processed AnnData object
adata = pipeline.adata


# -------------------
# ANNOTATION
# -------------------

annotator = CellAnnotator()

adata = annotator.annotate(adata)

cluster_map = annotator.cluster_annotations(
    adata
)


# -------------------
# MARKERS
# -------------------

marker_analyzer = MarkerAnalyzer()

adata = marker_analyzer.find_markers(
    adata
)

summary = marker_analyzer.summarize_clusters(
    adata,
    cluster_map
)


# -------------------
# PATHWAYS
# -------------------

pathway = PathwayAnalyzer()

for cluster in summary:

    markers = summary[cluster]["markers"]

    try:

        pathways = pathway.top_pathways(
            markers
        )

        summary[cluster]["pathways"] = pathways

    except Exception:

        summary[cluster]["pathways"] = None


# -------------------
# DISEASE INTERPRETATION
# -------------------

disease = DiseaseInterpreter()

for cluster in summary:

    pathways = summary[cluster]["pathways"]

    if pathways is not None:

        summary[cluster][
            "disease_interpretation"
        ] = disease.interpret_pathways(
            pathways["Term"]
            .head(5)
            .tolist()
        )

    else:

        summary[cluster][
            "disease_interpretation"
        ] = [
            "No pathway enrichment available."
        ]


# -------------------
# CONCLUSIONS
# -------------------

engine = ConclusionEngine()

for cluster in summary:

    pathways = ""

    if summary[cluster]["pathways"] is not None:

        pathways = " ".join(
            summary[cluster]["pathways"]["Term"]
            .head(5)
            .tolist()
        )

    summary[cluster][
        "conclusion"
    ] = engine.generate_conclusion(
        summary[cluster]["cell_type"],
        summary[cluster]["markers"],
        pathways
    )


# -------------------
# COMMUNICATION
# -------------------

communicator = CellCommunication()

communication = communicator.analyze(
    summary
)


# -------------------
# RESULTS
# -------------------

for cluster in summary:

    print("\n")
    print("=" * 60)

    print(
        f"Cluster {cluster}"
    )

    print(
        f"Cell Type: "
        f"{summary[cluster]['cell_type']}"
    )

    print(
        f"Markers: "
        f"{', '.join(summary[cluster]['markers'][:5])}"
    )

    print("\nDisease Interpretation:")

    for finding in summary[cluster][
        "disease_interpretation"
    ]:

        print(
            f"- {finding}"
        )

    print("\nConclusion:")

    print(
        summary[cluster]["conclusion"]
    )

    print("\nCommunication:")

    for signal in communication[
        cluster
    ]["signals"]:

        print(
            f"- {signal}"
        )
# -------------------
# FINAL REPORT
# -------------------

reporter = FinalReport()

report = reporter.generate(
    summary,
    communication
)

reporter.save(report)
