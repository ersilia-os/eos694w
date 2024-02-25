import sys
import os
import time
import csv
import json

import click
import torch
import numpy as np
from reinvent.config_parse import read_smiles_csv_file
from reinvent.runmodes.dtos import ChemistryHelpers
from reinvent.chemistry import Conversions
from reinvent.chemistry.library_design import BondMaker, AttachmentPoints
from reinvent.models.model_factory.sample_batch import SampleBatch
from reinvent.runmodes.samplers.run_sampling import filter_valid

from setup_sampler import setup_sampler

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
is_debug = sys.argv[3] == "True" if len(sys.argv) > 3 else False
log_file = output_file + ".json"


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


def filter_out_duplicate_molecules(sampled: SampleBatch) -> SampleBatch:
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


# my model
def my_model():
    batch_size = 100
    chemistry = ChemistryHelpers(Conversions(), BondMaker(), AttachmentPoints())
    sampler = setup_sampler(batch_size=batch_size, chemistry=chemistry)

    input_smiles = None
    num_input_smiles = 0

    if os.path.exists(input_file):
        input_smiles = read_smiles_csv_file(input_file, columns=0)
        num_input_smiles = len(input_smiles)
    else:
        click.echo(click.style(f"[INPUT_FILE]: {input_file} doesn't exist.", fg="red"))
        return

    if not os.path.exists(os.path.dirname(os.path.abspath(output_file))):
        click.echo(
            click.style(
                f"[OUTPUT_DIR]: {os.path.dirname(output_file)} doesn't exist.", fg="red"
            )
        )
        return

    start_time = time.time()

    if is_debug:
        click.echo(
            click.style(
                f"Starting sampling at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}",
                fg="green",
            )
        )
        click.echo(f"Total input smiles: {num_input_smiles}")

    with torch.no_grad():
        sampled = sampler.sample(input_smiles)
    end_time = time.time()

    if is_debug:
        click.echo(
            click.style(
                f"Time taken in seconds: {int(end_time - start_time)}", fg="green"
            )
        )

    # compute Tanimoto similarity between generated compounds and input compounds; return largest
    _, valid_idxs = chemistry.conversions.smiles_to_mols_and_indices(sampled.items2)
    valid_scores = sampler.calculate_tanimoto(input_smiles, sampled.items2)
    scores = [None] * len(sampled.items2)
    for i, j in enumerate(valid_idxs):
        scores[j] = valid_scores[i]

    sampled, scores = sort_sampled_molecules(sampled, scores)
    sampled = filter_valid(sampled)
    sampled = filter_out_duplicate_molecules(sampled)

    total_smiles = len(sampled.smilies)

    if is_debug:
        expected_num_smiles = batch_size * num_input_smiles
        click.echo(
            click.style(
                f"Total unique smiles generated: {total_smiles}, Expected: {expected_num_smiles}, Loss: {expected_num_smiles - total_smiles}"
            )
        )

    HEADER = [f"Smiles {x}" for x in range(total_smiles)]

    with open(output_file, "w", newline="") as fp:
        csv_writer = csv.writer(fp)
        # First Row: Header
        # Second Row: Generated Smiles (Output)
        # Third Row: Input Smiles
        csv_writer.writerows([HEADER, sampled.smilies, sampled.items1])

    with open(os.path.abspath(log_file), "w", newline="\n") as fp:
        log = {
            "start": start_time,
            "end": end_time,
            "input_smiles": num_input_smiles,
            "total": total_smiles,
            "expected": expected_num_smiles,
        }

        json.dump(log, fp)


# run model
my_model()
