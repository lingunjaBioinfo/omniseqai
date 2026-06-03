class FinalReport:

    def generate(
        self,
        cluster_summary,
        communication,
        de_findings=None
    ):

        report = []

        report.append(
            "========== OMNISEQAI REPORT ==========\n"
        )

        # --------------------------
        # CLUSTERS
        # --------------------------

        for cluster in cluster_summary:

            report.append(
                f"\nCluster {cluster}"
            )

            report.append(
                f"Cell Type: "
                f"{cluster_summary[cluster]['cell_type']}"
            )

            report.append(
                f"Markers: "
                f"{', '.join(cluster_summary[cluster]['markers'][:5])}"
            )

            report.append(
                f"Conclusion: "
                f"{cluster_summary[cluster]['conclusion']}"
            )

            if (
                "disease_interpretation"
                in cluster_summary[cluster]
            ):

                report.append(
                    "Disease Interpretation:"
                )

                for finding in (
                    cluster_summary[cluster]
                    ["disease_interpretation"]
                ):

                    report.append(
                        f"- {finding}"
                    )

            report.append(
                "Communication:"
            )

            for signal in communication[
                cluster
            ]["signals"]:

                report.append(
                    f"- {signal}"
                )

        # --------------------------
        # CONDITION ANALYSIS
        # --------------------------

        if de_findings is not None:

            report.append(
                "\n\n===== CONDITION ANALYSIS ====="
            )

            for finding in de_findings:

                report.append(
                    f"- {finding}"
                )

        return "\n".join(report)

    def save(
        self,
        report_text,
        filename="reports/final_analysis_report.txt"
    ):

        with open(
            filename,
            "w"
        ) as f:

            f.write(report_text)

        print(
            f"\nFinal report saved: {filename}"
        )
