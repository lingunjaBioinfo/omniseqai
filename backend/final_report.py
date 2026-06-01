class FinalReport:

    def generate(
        self,
        summary,
        communication
    ):

        report = []

        report.append(
            "OMNISEQAI ANALYSIS REPORT\n"
        )

        report.append(
            "=" * 60 + "\n"
        )

        for cluster in summary:

            report.append(
                f"\nCluster {cluster}\n"
            )

            report.append(
                f"Cell Type: "
                f"{summary[cluster]['cell_type']}\n"
            )

            report.append(
                "Markers:\n"
            )

            report.append(
                ", ".join(
                    summary[cluster]["markers"][:10]
                ) + "\n"
            )

            if (
                "disease_interpretation"
                in summary[cluster]
            ):

                report.append(
                    "\nDisease Interpretation:\n"
                )

                for item in summary[cluster][
                    "disease_interpretation"
                ]:

                    report.append(
                        f"- {item}\n"
                    )

            report.append(
                "\nConclusion:\n"
            )

            report.append(
                summary[cluster][
                    "conclusion"
                ] + "\n"
            )

            report.append(
                "\nCommunication:\n"
            )

            for signal in communication[
                cluster
            ]["signals"]:

                report.append(
                    f"- {signal}\n"
                )

            report.append(
                "\n" + "-" * 60 + "\n"
            )

        return "".join(report)

    def save(
        self,
        report,
        filename="final_analysis_report.txt"
    ):

        path = f"reports/{filename}"

        with open(
            path,
            "w"
        ) as f:

            f.write(report)

        print(
            f"\nFinal report saved: {path}"
        )
