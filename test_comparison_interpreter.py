import pandas as pd

from backend.comparison_interpreter import (
    ComparisonInterpreter
)

table = pd.DataFrame(
    {
        "Healthy": [0.33, 0.17, 0.50],
        "Disease": [0.50, 0.33, 0.17]
    },
    index=[
        "Monocytes",
        "NK cells",
        "T cells"
    ]
)

interpreter = ComparisonInterpreter()

results = interpreter.interpret(
    table
)

for item in results:

    print(item)
