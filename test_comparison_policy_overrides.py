import pandas as pd

from backend.comparison_policy import ComparisonPolicy


policy = ComparisonPolicy()

series = pd.Series(["ArmA", "ArmB", "Baseline", "ArmA", "ArmB", "Baseline"])

print("\n=== Automatic ===")
auto = policy.build_plan(series)
print("All groups:", auto.all_groups)
print("Baseline:", auto.baseline)
print("Comparisons:", auto.comparisons)

print("\n=== Explicit control/case ===")
explicit = policy.build_plan(series, control="Baseline", case="ArmA")
print("All groups:", explicit.all_groups)
print("Baseline:", explicit.baseline)
print("Comparisons:", explicit.comparisons)

print("\n=== Explicit pairs ===")
paired = policy.build_plan(series, pairs=[("Baseline", "ArmA"), ("Baseline", "ArmB")])
print("All groups:", paired.all_groups)
print("Baseline:", paired.baseline)
print("Comparisons:", paired.comparisons)
