import os
import time

import click
import torch
from reinvent.runmodes import samplers, create_adapter
from reinvent.runmodes.dtos import ChemistryHelpers
from reinvent.chemistry import TransformationTokens, Conversions
from reinvent.chemistry.library_design import BondMaker, AttachmentPoints
from reinvent.runmodes.samplers.run_sampling import filter_valid


from utils import (
    sort_sampled_molecules,
    filter_out_duplicate_molecules,
    pad_smiles,
    make_list_into_lists_of_n,
)


class Mol2MolSimilaritySampler:
    ROOT = os.path.dirname(os.path.abspath(__file__))
    CHECKPOINT = os.path.join(ROOT, "..", "..", "checkpoints")
    MODEL = os.path.realpath(os.path.join(CHECKPOINT, "mol2mol_similarity.prior"))

    def __init__(self, batch_size: int):
        self.batch_size = batch_size
        self.chemistry = ChemistryHelpers(
            Conversions(), BondMaker(), AttachmentPoints()
        )

        # Constants
        temperature = 1.0
        unique_sequences = True
        isomeric = True
        sample_strategy = "beamsearch"
        tokens = TransformationTokens()
        randomize_smiles = (
            False  # For Mol2Mol reinvent 4 sets randomize_smiles to False.
        )

        # Creating adapter
        agent, _, _ = create_adapter(self.MODEL, "inference", torch.device("cpu"))

        # Creating sampler.
        self.sampler = samplers.Mol2MolSampler(
            agent,
            batch_size=self.batch_size,
            sample_strategy=sample_strategy,
            isomeric=isomeric,
            randomize_smiles=randomize_smiles,
            unique_sequences=unique_sequences,
            chemistry=self.chemistry,
            temperature=temperature,
            tokens=tokens,
        )

    def generate(
        self, input_smiles: "list[str]", is_debug: bool
    ) -> "tuple[list[list[str]], list[str], dict]":
        start_time = time.time()
        num_input_smiles = len(input_smiles)

        if is_debug:
            click.echo(
                click.style(
                    f"Starting sampling at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}",
                    fg="green",
                )
            )
            click.echo(f"Total input smiles: {num_input_smiles}")

        with torch.no_grad():
            sampled = self.sampler.sample(input_smiles)
            end_time = time.time()

            if is_debug:
                click.echo(
                    click.style(
                        f"Time taken in seconds: {int(end_time - start_time)}",
                        fg="green",
                    )
                )

        # compute Tanimoto similarity between generated compounds and input compounds; return largest
        _, valid_idxs = self.chemistry.conversions.smiles_to_mols_and_indices(
            sampled.items2
        )
        valid_scores = self.sampler.calculate_tanimoto(input_smiles, sampled.items2)
        scores = [-1] * len(sampled.items2)
        for i, j in enumerate(valid_idxs):
            scores[j] = valid_scores[i]

        sampled, scores = sort_sampled_molecules(sampled, scores)

        sampled = filter_valid(sampled)
        sampled = filter_out_duplicate_molecules(sampled, is_debug=is_debug)

        total_smiles = len(sampled.smilies)

        expected_num_smiles = self.batch_size * num_input_smiles

        if is_debug:
            click.echo(
                click.style(
                    f"Total unique smiles generated: {total_smiles}, Expected: {expected_num_smiles}, Loss: {expected_num_smiles - total_smiles}"
                )
            )

        flatten_outputs = pad_smiles(
            sampled, input_smiles=input_smiles, target_length=self.batch_size
        )

        output_smiles = make_list_into_lists_of_n(flatten_outputs, num_input_smiles)

        log = {
            "start": start_time,
            "end": end_time,
            "input_smiles": input_smiles,
            "total": total_smiles,
            "expected": expected_num_smiles,
        }

        return (output_smiles, flatten_outputs, log)
