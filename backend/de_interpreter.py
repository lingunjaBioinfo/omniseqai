class DEInterpreter:

    def interpret(
        self,
        de_results
    ):

        genes = (
            de_results["names"]
            .astype(str)
            .str.upper()
            .tolist()
        )

        findings = []

        # --------------------------
        # T CELLS
        # --------------------------

        tcell_markers = {
            "IL32",
            "LTB",
            "CD3D",
            "CD3E",
            "TRBC1",
            "TRBC2",
            "LTB"
        }

        if any(g in genes for g in tcell_markers):

            findings.append(
                "Enhanced T-cell activation detected."
            )

        # --------------------------
        # NK CELLS
        # --------------------------

        nk_markers = {
            "NKG7",
            "GNLY",
            "PRF1",
            "GZMB"
        }

        if any(g in genes for g in nk_markers):

            findings.append(
                "Increased cytotoxic immune activity."
            )

        # --------------------------
        # MONOCYTES
        # --------------------------

        mono_markers = {
            "LYZ",
            "S100A8",
            "S100A9",
            "FCN1"
        }

        if any(g in genes for g in mono_markers):

            findings.append(
                "Monocyte-driven inflammatory response."
            )

        # --------------------------
        # INTERFERON
        # --------------------------

        interferon_markers = {
            "ISG15",
            "IFIT1",
            "IFIT2",
            "IFIT3",
            "MX1",
            "OAS1"
        }

        if any(g in genes for g in interferon_markers):

            findings.append(
                "Interferon signaling suggests antiviral activity."
            )

        # --------------------------
        # DEFAULT
        # --------------------------

        if len(findings) == 0:

            findings.append(
                "No strong biological pattern detected."
            )

        return findings
