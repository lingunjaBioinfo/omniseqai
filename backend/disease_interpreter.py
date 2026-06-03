class DiseaseInterpreter:

    def interpret(self, markers):
        markers = set(str(m).upper() for m in markers)
        findings = []

        interferon = {
            "IFI27", "IFITM3", "IFI6", "ISG15", "MX1", "OAS1", "OASL", "IFI44L"
        }
        inflammatory = {
            "S100A8", "S100A9", "FCN1", "CXCL8", "IL1B", "NCF1"
        }
        cytotoxic = {
            "NKG7", "GNLY", "PRF1", "GZMB", "CTSW"
        }
        bcell = {
            "MS4A1", "CD79A", "CD79B", "CD74"
        }

        if len(interferon & markers) >= 2:
            findings.append("Type-I interferon response.")

        if len(inflammatory & markers) >= 2:
            findings.append("Inflammatory monocyte activation.")

        if len(cytotoxic & markers) >= 2:
            findings.append("Cytotoxic immune activation.")

        if len(bcell & markers) >= 2:
            findings.append("Activated B-cell response.")

        mt_count = sum(g.startswith("MT-") for g in markers)
        if mt_count >= 3:
            findings.append("Low-quality or stressed cell population.")

        if not findings:
            findings.append("No strong disease signature identified.")

        return findings
