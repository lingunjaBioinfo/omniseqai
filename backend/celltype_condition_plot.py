import matplotlib.pyplot as plt


class CellTypeConditionPlot:

    def create(
        self,
        proportions,
        output_file="outputs/celltype_conditions.png"
    ):

        proportions.plot(
            kind="bar"
        )

        plt.ylabel("Proportion")

        plt.title(
            "Cell Type Composition"
        )

        plt.tight_layout()

        plt.savefig(output_file)

        plt.close()

        print(
            f"\nPlot saved: "
            f"{output_file}"
        )
