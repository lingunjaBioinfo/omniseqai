import pandas as pd

from backend.comparison_policy import ComparisonPolicy


policy = ComparisonPolicy()

examples = {
    "COVID dataset": ["COVID", "Healthy"],
    "Vaccination dataset": ["Vaccinated", "Healthy", "COVID"],
    "Tumor dataset": ["Tumor", "Normal"],
    "Treatment dataset": ["DrugA", "DrugB", "Control"],
    "No obvious control": ["Group1", "Group2", "Group3"],
}

for name, labels in examples.items():
    series = pd.Series(labels)

    plan = policy.build_plan(series)

    print("\n" + "=" * 70)
    print(name)
    print("=" * 70)
    print("All groups:", plan.all_groups)
    print("Baseline:", plan.baseline)
    print("Control-like groups:", plan.control_like_groups)
    print("Comparisons:", plan.comparisons)
