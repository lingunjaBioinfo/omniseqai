import matplotlib.pyplot as plt
import scanpy as sc


class ConditionVisualization:

    def celltype_proportions(
        self,
        proportions,
        output_file="outputs/celltype_proportions.png"
    ):

        proportions.plot.bar()

        plt.title(
            "Cell Type Proportions by Condition"
        )

        plt.ylabel("Proportion")

        plt.tight_layout()

        plt.savefig(output_file)

        plt.close()

        print(
            f"Saved: {output_file}"
        )

    def condition_umap(
        self,
        adata,
        condition_key="condition",
        output_file="outputs/condition_umap.png"
    ):

        sc.pl.umap(
            adata,
            color=condition_key,
            show=False
        )

        plt.savefig(output_file)

        plt.close()

        print(
            f"Saved: {output_file}"
        )
