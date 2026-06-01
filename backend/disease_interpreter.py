class DiseaseInterpreter:

    def interpret_pathways(
        self,
        pathways
    ):

        pathway_text = " ".join(
            pathways
        ).lower()

        findings = []

        score = {
            "viral": 0,
            "inflammation": 0,
            "adaptive": 0
        }

        if "interferon" in pathway_text:
            score["viral"] += 2

        if "viral" in pathway_text:
            score["viral"] += 2

        if "cytokine" in pathway_text:
            score["inflammation"] += 1

        if "inflammatory" in pathway_text:
            score["inflammation"] += 2

        if "t cell activation" in pathway_text:
            score["adaptive"] += 2

        if "alpha-beta t cell activation" in pathway_text:
            score["adaptive"] += 1

        if score["viral"] >= 2:

            findings.append(
                "Possible antiviral immune response."
            )

        if score["inflammation"] >= 2:

            findings.append(
                "Inflammatory activation detected."
            )

        if score["adaptive"] >= 2:

            findings.append(
                "Activated adaptive immune response detected."
            )

        if len(findings) == 0:

            findings.append(
                "No strong disease signature identified."
            )

        return findings
