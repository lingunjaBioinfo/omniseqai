class ComparisonInterpreter:

    def interpret(
        self,
        proportions
    ):

        conditions = list(
            proportions.columns
        )

        if len(conditions) != 2:
            return [
                "Exactly two conditions required."
            ]

        control = conditions[0]
        disease = conditions[1]

        findings = []

        for cell_type in proportions.index:

            control_value = (
                proportions.loc[
                    cell_type,
                    control
                ]
            )

            disease_value = (
                proportions.loc[
                    cell_type,
                    disease
                ]
            )

            if disease_value > control_value:

                findings.append(
                    f"{cell_type} increased "
                    f"from {control_value:.2f} "
                    f"to {disease_value:.2f}"
                )

            elif disease_value < control_value:

                findings.append(
                    f"{cell_type} decreased "
                    f"from {control_value:.2f} "
                    f"to {disease_value:.2f}"
                )

        return findings
