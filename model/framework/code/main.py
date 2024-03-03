import sys
import os
import time
import csv
import json

import click
import torch
from reinvent.config_parse import read_smiles_csv_file
from reinvent.runmodes.dtos import ChemistryHelpers
from reinvent.chemistry import Conversions
from reinvent.chemistry.library_design import BondMaker, AttachmentPoints
from reinvent.runmodes.samplers.run_sampling import filter_valid

from setup_sampler import setup_sampler
from utils import sort_sampled_molecules, filter_out_duplicate_molecules, pad_smiles

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
is_debug = sys.argv[3] == "True" if len(sys.argv) > 3 else False
log_file = output_file + ".json"


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
    scores = [-1] * len(sampled.items2)
    for i, j in enumerate(valid_idxs):
        scores[j] = valid_scores[i]

    sampled, scores = sort_sampled_molecules(sampled, scores)
    sampled = filter_valid(sampled)
    sampled = filter_out_duplicate_molecules(sampled, is_debug=is_debug)

    total_smiles = len(sampled.smilies)

    if is_debug:
        expected_num_smiles = batch_size * num_input_smiles
        click.echo(
            click.style(
                f"Total unique smiles generated: {total_smiles}, Expected: {expected_num_smiles}, Loss: {expected_num_smiles - total_smiles}"
            )
        )

    output_smiles = pad_smiles(
        sampled, input_smiles=input_smiles, target_length=batch_size
    )
    HEADER = [f"smi_{x}" for x in range(num_input_smiles * batch_size)]

    with open(output_file, "w", newline="") as fp:
        csv_writer = csv.writer(fp)
        # First Row: Header
        # Second Row: Generated Smiles (Output)
        csv_writer.writerows([HEADER, output_smiles])

    with open(os.path.abspath(log_file), "w", newline="\n") as fp:
        log = {
            "start": start_time,
            "end": end_time,
            "input_smiles": input_smiles,
            "total": total_smiles,
            "expected": expected_num_smiles,
        }

        json.dump(log, fp)


# run model
my_model()
