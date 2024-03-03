import click
import numpy as np
from reinvent.models.model_factory.sample_batch import SampleBatch
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem


def are_smiles_same(smiles1, smiles2):
    """To check whether given smiles are same."""
    # Parse SMILES strings to obtain RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Check if both SMILES strings are valid
    if mol1 is None or mol2 is None:
        return False  # Unable to parse one or both SMILES strings

    # Generate Morgan fingerprints for each molecule
    fp1 = Chem.AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)

    # Compare fingerprints to assess similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

    # Determine whether the molecules are considered the same
    return similarity == 1.0


def sort_sampled_molecules(
    sampled: SampleBatch, scores: "list[float]"
) -> "tuple[SampleBatch, list[float]]":
    """Sort sampled molecules using tanimoto in descending order."""
    items1 = np.array(sampled.items1)
    items2 = np.array(sampled.items2)
    smilies = np.array(sampled.smilies)
    states = np.array(sampled.states)
    _scores = np.array(scores)

    sorted_indices = np.argsort(_scores)[::-1]

    sampled_batch = SampleBatch(
        items1=list(items1[sorted_indices]),
        items2=list(items2[sorted_indices]),
        nlls=sampled.nlls[sorted_indices.copy()],
        smilies=list(smilies[sorted_indices]),
        states=list(states[sorted_indices]),
    )

    return sampled_batch, list(_scores[sorted_indices])


def filter_out_duplicate_molecules(sampled: SampleBatch, is_debug: bool) -> SampleBatch:
    """Filter out duplicate molecules from the sampled molecules.
    It also remove the output molecules if it is similar to input molecules.

    `sampled.items1` contains input smiles.
    `sampled.smilies` contains output smiles.
    """

    seen = {}
    items1 = []
    items2 = []
    states = []
    smilies = []
    nll_indices = []

    for item1, item2, smile, nll_index, state in zip(
        sampled.items1,
        sampled.items2,
        sampled.smilies,
        range(len(sampled.items1)),
        sampled.states,
    ):
        seen[item1] = item1

        if smile in seen:
            if is_debug:
                click.echo(
                    click.style(
                        f"Removing {smile}, as it is a duplicate entry.", fg="yellow"
                    )
                )
            continue

        seen[smile] = smile

        items1.append(item1)
        items2.append(item2)
        smilies.append(smile)
        states.append(state)
        nll_indices.append(nll_index)

    return SampleBatch(
        items1=items1,
        items2=items2,
        states=states,
        smilies=smilies,
        nlls=sampled.nlls[nll_indices],
    )


def pad_smiles(
    sampled: SampleBatch, input_smiles: "list[str]", target_length=100
) -> "list[str]":
    """For a given input smiles, it is not possible to get the expected
    number (target_length) of output smiles (samples). This will cause
    problems in the downstream process. To mitigate this, the function
    pads the output with empty strings to match the target_length.
    """

    output_smiles: dict[list[str]] = {}

    # Initialize output_smiles with empty lists for each input smile
    output_smiles = {smile: [] for smile in input_smiles}

    # Smiles with similar tanimoto.
    similar_smiles = {}

    # Populate output_smiles with sampled smiles
    for idx, seq in enumerate(sampled.items1):
        if seq in output_smiles:
            output_smiles[seq].append(sampled.smilies[idx])
        else:
            # REINVENT4 do some modification in the input smiles. Like if the
            # input smile is `CC(=O)Oc1ccccc1C(O)=O` then it will convert it
            # to `CC(=O)Oc1ccccc1C(=O)O`. Both of them are same, however both
            # strings are not same. This else condition will take care of this
            # edge case.
            if seq in similar_smiles:
                output_smiles[similar_smiles[seq]].append(sampled.smilies[idx])
            else:
                # Try to find similar smiles in input_smiles.
                for smile in input_smiles:
                    if are_smiles_same(smile, seq):
                        output_smiles[smile].append(sampled.smilies[idx])
                        similar_smiles[seq] = smile
                        break

    output = []

    # Construct the output list
    for smile in input_smiles:
        output_smile = output_smiles[smile]
        output_smile_length = len(output_smile)
        if output_smile_length == target_length:
            output.extend(output_smile)
        else:
            output.extend(output_smile)
            padding = [""] * (target_length - output_smile_length)
            output.extend(padding)

    return output
