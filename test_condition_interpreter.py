from backend.condition_interpreter import (
    ConditionInterpreter
)

terms = [
    "T Cell Activation",
    "Interferon Gamma Signaling"
]

interpreter = ConditionInterpreter()

results = interpreter.interpret(
    terms
)

for item in results:
    print(item)
