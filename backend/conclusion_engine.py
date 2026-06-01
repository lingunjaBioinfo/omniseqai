class ConclusionEngine:

    def generate_conclusion(
        self,
        cell_type,
        markers,
        pathways
    ):

        text = []

        text.append(
            f"The dominant population "
            f"consists of {cell_type}."
        )

        if any(
            x in pathways.lower()
            for x in [
                "t cell activation",
                "adaptive immune"
            ]
        ):
            text.append(
                "Evidence suggests activation "
                "of adaptive immunity."
            )

        if any(
            x in pathways.lower()
            for x in [
                "interferon",
                "viral"
            ]
        ):
            text.append(
                "This pattern may indicate "
                "a response to viral infection."
            )

        if any(
            x in pathways.lower()
            for x in [
                "inflammatory"
            ]
        ):
            text.append(
                "Inflammatory signaling appears "
                "to be elevated."
            )

        return " ".join(text)
