class ConclusionEngine:

    def generate(
        self,
        cell_type
    ):

        cell_type_lower = (
            cell_type.lower()
        )

        if "monocyte" in cell_type_lower:

            return (
                f"{cell_type} population "
                f"driving innate immune responses."
            )

        elif "nk" in cell_type_lower:

            return (
                f"{cell_type} population "
                f"consistent with antiviral "
                f"cytotoxic activity."
            )

        elif "t cell" in cell_type_lower:

            return (
                f"{cell_type} population "
                f"involved in adaptive immunity."
            )

        elif "b cell" in cell_type_lower:

            return (
                f"{cell_type} population "
                f"associated with humoral immunity."
            )

        elif "dendritic" in cell_type_lower:

            return (
                f"{cell_type} population "
                f"likely involved in antigen "
                f"presentation."
            )

        return (
            f"Dominant population: "
            f"{cell_type}"
        )
