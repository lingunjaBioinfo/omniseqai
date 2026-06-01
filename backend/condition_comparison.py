import pandas as pd


class ConditionComparison:

    def compare_celltypes(
        self,
        adata,
        condition_column="condition",
        celltype_column="cell_type"
    ):

        proportions = pd.crosstab(
            adata.obs[celltype_column],
            adata.obs[condition_column],
            normalize="columns"
        )

        return proportions

    def summarize_changes(
        self,
        proportions
    ):

        print(
            "\n===== CELL TYPE PROPORTIONS =====\n"
        )

        print(
            proportions
        )

        return proportions
