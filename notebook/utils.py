import rdkit
import matplotlib.pyplot as plt
import csv
import os
import numpy as np
import json
import unicodedata
import string

valid_filename_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
char_limit = 255


def molecule_from_smiles_img(smile):
    return rdkit.Chem.Draw.MolToImage(rdkit.Chem.MolFromSmiles(smile))


def read_molecule(input_path):
    with open(input_path) as fp:
        smi = fp.read()
        return smi


def clean_filename(filename, whitelist=valid_filename_chars, replace="-"):
    # replace spaces
    for r in replace:
        filename = filename.replace(r, "_")

    # keep only valid ascii chars
    cleaned_filename = (
        unicodedata.normalize("NFKD", filename).encode("ASCII", "ignore").decode()
    )

    # keep only whitelisted chars
    cleaned_filename = "".join(c for c in cleaned_filename if c in whitelist)

    return cleaned_filename[:char_limit]


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

    inputs = []

    with open(f"{output_smiles_path}.json") as logfile:
        log = json.load(logfile)
        batch_length = log["expected"] // len(log["input_smiles"])

        for smile in log["input_smiles"]:
            smiles = [smile] * batch_length
            inputs.extend(smiles)

    inputs = np.array(inputs)
    outputs = np.array(rows[0])

    perm = np.random.permutation(len(inputs))

    image_output_dir = os.path.dirname(output_smiles_path)

    for input_smile, output_smile in zip(inputs[perm], outputs[perm]):
        if not output_smile:
            continue

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
            os.path.join(image_output_dir, clean_filename(output_smile)),
            dpi=300,
            bbox_inches="tight",
        )

        plt.show()
