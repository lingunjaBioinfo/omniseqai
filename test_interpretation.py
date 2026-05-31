from backend.interpretation import (
    BiologicalInterpreter
)

interpreter = BiologicalInterpreter()

summary = interpreter.generate_summary(
    cluster=1,
    cell_type="Classical monocytes",
    markers=[
        "LYZ",
        "S100A8",
        "S100A9",
        "FCN1",
        "TYROBP"
    ]
)

print(summary)
