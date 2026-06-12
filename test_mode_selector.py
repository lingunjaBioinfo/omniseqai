from backend.mode_selector import ModeSelector

selector = ModeSelector()

profile_condition = {
    "condition_column": "disease",
    "cell_type_column": "cell_type",
    "sample_column": "sample_id",
    "patient_column": "donor_id",
}

profile_exploratory = {
    "condition_column": None,
    "cell_type_column": "cell_type",
    "sample_column": None,
    "patient_column": None,
}

print(selector.choose(profile_condition, mode="auto"))
print(selector.choose(profile_exploratory, mode="auto"))
print(selector.choose(profile_condition, mode="exploratory"))
print(selector.choose(profile_condition, mode="condition"))
