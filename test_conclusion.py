from backend.conclusion_engine import ConclusionEngine

engine = ConclusionEngine()

result = engine.generate_conclusion(
    "T cells",
    ["CD3D", "IL7R"],
    """
    T Cell Activation
    Adaptive Immune Response
    """
)

print(result)
