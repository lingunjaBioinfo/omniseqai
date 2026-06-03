import os

import numpy as np
import matplotlib.pyplot as plt


class VolcanoPlot:

    def create(
        self,
        de_results,
        output_file="outputs/volcano_plot.png",
        title="Condition Differential Expression Volcano Plot"
    ):

        if de_results is None or de_results.empty:
            print("No DE results available for volcano plot.")
            return

        df = de_results.copy()

        df["logfoldchanges"] = np.asarray(
            df["logfoldchanges"],
            dtype=float
        )
        df["pvals_adj"] = np.asarray(
            df["pvals_adj"],
            dtype=float
        )

        df["minus_log10_p"] = -np.log10(
            np.clip(df["pvals_adj"], 1e-300, 1.0)
        )

        sig = (
            (df["pvals_adj"] < 0.05) &
            (df["logfoldchanges"].abs() > 1.0)
        )

        plt.figure(figsize=(9, 7))

        plt.scatter(
            df.loc[~sig, "logfoldchanges"],
            df.loc[~sig, "minus_log10_p"],
            s=10,
            alpha=0.5,
            label="Not significant"
        )

        plt.scatter(
            df.loc[sig, "logfoldchanges"],
            df.loc[sig, "minus_log10_p"],
            s=12,
            alpha=0.8,
            label="Significant"
        )

        plt.axvline(1.0, linestyle="--")
        plt.axvline(-1.0, linestyle="--")
        plt.axhline(-np.log10(0.05), linestyle="--")

        plt.xlabel("Log fold change")
        plt.ylabel("-log10(adjusted p-value)")
        plt.title(title)
        plt.legend(frameon=False)
        plt.tight_layout()

        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"\nVolcano plot saved: {output_file}")
