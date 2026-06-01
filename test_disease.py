from backend.disease_interpreter import DiseaseInterpreter


pathways = [
    "T Cell Activation",
    "Response To Interleukin-4",
    "Alpha-Beta T Cell Activation"
]

interpreter = DiseaseInterpreter()

results = interpreter.interpret_pathways(
    pathways
)

print(
    "\n===== DISEASE INTERPRETATION =====\n"
)

for item in results:

    print(
        f"- {item}"
    )
