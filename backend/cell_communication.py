class CellCommunication:

    def __init__(self):

        self.communication_rules = {

            "T cells": [
                "Adaptive immune signaling",
                "T-cell activation"
            ],

            "B cells": [
                "Antibody-mediated signaling",
                "Antigen presentation"
            ],

            "Classical monocytes": [
                "Inflammatory cytokine signaling",
                "Innate immune activation"
            ],

            "Non-classical monocytes": [
                "Immune surveillance",
                "Chemokine signaling"
            ],

            "CD16+ NK cells": [
                "Cytotoxic signaling",
                "Immune cell killing"
            ],

            "NK cells": [
                "Cytotoxic signaling"
            ],

            "Megakaryocytes/platelets": [
                "Platelet activation",
                "Coagulation signaling"
            ]
        }

    def analyze(
        self,
        cluster_summaries
    ):

        results = {}

        for cluster, info in cluster_summaries.items():

            cell_type = info["cell_type"]

            signals = self.communication_rules.get(
                cell_type,
                ["Unknown signaling"]
            )

            results[cluster] = {

                "cell_type": cell_type,
                "signals": signals
            }

        return results
