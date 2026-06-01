from backend.cell_communication import CellCommunication


sample_clusters = {

    0: {
        "cell_type": "T cells"
    },

    1: {
        "cell_type": "Classical monocytes"
    },

    2: {
        "cell_type": "CD16+ NK cells"
    }
}


communicator = CellCommunication()

results = communicator.analyze(
    sample_clusters
)

print(
    "\n===== CELL COMMUNICATION =====\n"
)

for cluster, info in results.items():

    print(
        f"Cluster {cluster}"
    )

    print(
        f"Cell Type: {info['cell_type']}"
    )

    print(
        "Potential Signals:"
    )

    for signal in info["signals"]:

        print(
            f"- {signal}"
        )

    print()
