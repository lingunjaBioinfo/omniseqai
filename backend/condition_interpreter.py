class ConditionInterpreter:

    IMMUNE_RULES = [
        (
            ["interferon", "antiviral", "virus", "viral"],
            "Possible antiviral interferon response."
        ),
        (
            ["cytokine", "chemokine", "inflammatory", "tnf", "il1", "nf-kappa"],
            "Inflammatory signaling increased."
        ),
        (
            ["t cell activation", "alpha-beta t cell", "adaptive immune"],
            "Enhanced adaptive immune activation."
        ),
        (
            ["antigen presentation", "mhc", "cd74", "hla-dra", "hla-dr"],
            "Antigen presentation activity detected."
        ),
        (
            ["nk cell", "cytotoxic", "granzyme", "perforin"],
            "Cytotoxic immune activity detected."
        ),
    ]

    HOUSEKEEPING_HINTS = [
        "translation",
        "ribosome",
        "spliceosome",
        "protein-containing complex assembly",
        "cellular respiration",
        "oxidative phosphorylation",
    ]

    def interpret(self, pathway_terms):
        """
        pathway_terms can be:
        - a list of strings
        - a pandas DataFrame with a 'Term' column
        """

        if pathway_terms is None:
            return ["No strong condition-specific pathway signature identified."]

        if hasattr(pathway_terms, "columns") and "Term" in pathway_terms.columns:
            terms = pathway_terms["Term"].astype(str).tolist()
        else:
            terms = [str(x) for x in pathway_terms]

        text = " ".join(terms).lower()
        findings = []

        for hints, label in self.IMMUNE_RULES:
            if any(h in text for h in hints):
                findings.append(label)

        if any(h in text for h in self.HOUSEKEEPING_HINTS):
            findings.append("Housekeeping or translational stress signature present.")

        # Remove duplicates while preserving order.
        seen = set()
        deduped = []
        for item in findings:
            if item not in seen:
                seen.add(item)
                deduped.append(item)

        if not deduped:
            deduped.append("No strong condition-specific pathway signature identified.")

        return deduped
