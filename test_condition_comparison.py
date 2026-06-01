import pandas as pd

from backend.condition_comparison import (
    ConditionComparison
)


data = {

    "cell_type": [

        "T cells",
        "T cells",
        "T cells",

        "Monocytes",
        "Monocytes",

        "NK cells",

        "T cells",

        "Monocytes",
        "Monocytes",
        "Monocytes",

        "NK cells",
        "NK cells"
    ],

    "condition": [

        "Healthy",
        "Healthy",
        "Healthy",

        "Healthy",
        "Healthy",

        "Healthy",

        "Disease",

        "Disease",
        "Disease",
        "Disease",

        "Disease",
        "Disease"
    ]
}

obs = pd.DataFrame(data)

class MockAdata:
    pass

adata = MockAdata()

adata.obs = obs

comparison = ConditionComparison()

table = comparison.compare_celltypes(
    adata
)

comparison.summarize_changes(
    table
)
