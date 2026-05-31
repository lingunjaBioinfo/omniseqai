class BiologicalInterpreter:

    def generate_summary(
        self,
        cluster,
        cell_type,
        markers
    ):

        top_markers = ", ".join(
            markers[:5]
        )

        return (
            f"Cluster {cluster} was annotated as "
            f"{cell_type}. "
            f"Key marker genes include "
            f"{top_markers}. "
            f"This transcriptional profile is "
            f"consistent with {cell_type}."
        )
