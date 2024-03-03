import rdkit
import matplotlib.pyplot as plt
import csv
import os
import numpy as np


def molecule_from_smiles_img(smile):
    return rdkit.Chem.Draw.MolToImage(rdkit.Chem.MolFromSmiles(smile))


def read_molecule(input_path):
    with open(input_path) as fp:
        smi = fp.read()
        smi = rdkit.Chem.Draw.MolToImage(rdkit.Chem.MolFromSmiles(smi))
        return smi


def visualize_random_output_molecules(output_smiles_path):
    rows = []
    with open(output_smiles_path, newline="") as csvfile:
        reader = csv.reader(csvfile)
        count = 0
        for row in reader:
            if count == 0:
                # Header
                count += 1
                continue
            rows.append(row)

    visited_input_smiles = {}

    inputs = np.array(rows[1])
    outputs = np.array(rows[0])

    perm = np.random.permutation(len(inputs))

    for input_smile, output_smile in zip(inputs[perm], outputs[perm]):
        current_input_smile = input_smile

        if current_input_smile in visited_input_smiles:
            if visited_input_smiles[current_input_smile] > 4:
                continue
            else:
                visited_input_smiles[current_input_smile] += 1
        else:
            visited_input_smiles[current_input_smile] = 1

        count += 1

        fig, axs = plt.subplots(
            1, 2, figsize=(10, 5)
        )  # Create a subplot with 1 row and 2 columns

        # Draw the input molecule
        axs[0].imshow(molecule_from_smiles_img(input_smile))
        axs[0].set_title("Input")
        axs[0].axis("off")  # Hide axes

        # Draw the output molecule
        axs[1].imshow(molecule_from_smiles_img(output_smile))
        axs[1].set_title("Output")
        axs[1].axis("off")  # Hide axes

        plt.savefig(
            os.path.join(output_smiles_path, "..", output_smile),
            dpi=300,
            bbox_inches="tight",
        )

        plt.show()