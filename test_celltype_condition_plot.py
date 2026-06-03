import pandas as pd

from backend.condition_comparison import (
    ConditionComparison
)

from backend.celltype_condition_plot import (
    CellTypeConditionPlot
)


obs = pd.DataFrame(
    {
        "cell_type": [
            "T cells",
            "T cells",
            "Monocytes",
            "NK cells",
            "Monocytes",
            "Monocytes"
        ],
        "condition": [
            "Healthy",
            "Healthy",
            "Healthy",
            "Disease",
            "Disease",
            "Disease"
        ]
    }
)

class Dummy:
    pass

adata = Dummy()

adata.obs = obs

comparison = ConditionComparison()

props = comparison.compare_celltypes(
    adata
)

plotter = CellTypeConditionPlot()

plotter.create(props)
