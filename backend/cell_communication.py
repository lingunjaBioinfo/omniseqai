class CellCommunication:

    def analyze(self, summary):
        communication = {}

        for cluster in summary:
            markers = set(str(x).upper() for x in summary[cluster]["markers"])
            signals = []

            if {"CXCL8", "CXCL2", "CCL2", "CCL3"} & markers:
                signals.append("Inflammatory chemokine signaling")

            if {"IFI27", "IFITM3", "ISG15", "MX1"} & markers:
                signals.append("Interferon signaling")

            if {"GNLY", "NKG7", "PRF1", "GZMB"} & markers:
                signals.append("Cytotoxic signaling")

            if {"HLA-DRA", "HLA-DRB1", "CD74"} & markers:
                signals.append("Antigen presentation")

            if not signals:
                signals.append("Unknown signaling")

            communication[cluster] = {"signals": signals}

        return communication
