import os
import sys
import time

import torch
from rdkit import Chem
from reinvent.runmodes import samplers, create_adapter
from reinvent.runmodes.dtos import ChemistryHelpers
from reinvent.chemistry import TransformationTokens, Conversions
from reinvent.chemistry.library_design import BondMaker, AttachmentPoints
from reinvent.runmodes.samplers.run_sampling import filter_valid

from utils import is_same_2d_structure

N_OUTPUTS = 100  # columns of run_columns.csv
DRAWS_PER_ROUND = 150  # sequences drawn by one call of the sampler
MAX_ROUNDS = 5  # give up on an input after this many rounds
TEMPERATURE = 1.0  # REINVENT 4's default for sampling


class Mol2MolSimilaritySampler:
    """Samples analogues of a molecule with a similarity-trained Mol2Mol prior of REINVENT 4.

    Each round draws DRAWS_PER_ROUND sequences by multinomial sampling (unseeded, so every call
    gives a different set). Rounds are repeated until N_OUTPUTS distinct molecules that are not the
    input (stereochemistry ignored) are collected, or MAX_ROUNDS is reached; the row is then
    padded with empty cells. Nothing is ranked: the molecules stay in draw order.
    """

    ROOT = os.path.dirname(os.path.abspath(__file__))
    CHECKPOINT = os.path.join(ROOT, "..", "..", "checkpoints")
    MODEL = os.path.realpath(os.path.join(CHECKPOINT, "mol2mol_similarity.prior"))

    def __init__(self):
        self.chemistry = ChemistryHelpers(
            Conversions(), BondMaker(), AttachmentPoints()
        )
        agent, _, _ = create_adapter(self.MODEL, "inference", torch.device("cpu"))
        self.sampler = samplers.Mol2MolSampler(
            agent,
            batch_size=DRAWS_PER_ROUND,
            sample_strategy="multinomial",
            isomeric=True,
            randomize_smiles=False,  # REINVENT 4 sets it to False for Mol2Mol
            unique_sequences=True,
            chemistry=self.chemistry,
            temperature=TEMPERATURE,
            tokens=TransformationTokens(),
        )

    def sample_one(self, smiles: str) -> "list[str]":
        """Up to N_OUTPUTS distinct analogues of `smiles`, in draw order."""
        pool = {}
        for _ in range(MAX_ROUNDS):
            with torch.no_grad():
                sampled = filter_valid(self.sampler.sample([smiles]))
            for smile in sampled.smilies:
                mol = Chem.MolFromSmiles(smile)
                if mol is None or is_same_2d_structure(smile, smiles):
                    continue
                # canonical isomeric SMILES as the key keeps one entry per molecule
                pool.setdefault(Chem.MolToSmiles(mol), None)
            if len(pool) >= N_OUTPUTS:
                break
        return list(pool)[:N_OUTPUTS]

    def generate(
        self, input_smiles: "list[str]", is_debug: bool
    ) -> "tuple[list[list[str]], list[str], dict]":
        start_time = time.time()
        output_smiles = []
        total_smiles = 0
        for smi in input_smiles:
            try:
                row = self.sample_one(smi)
            except Exception as exc:  # one failing molecule must not sink the batch
                print(f"Mol2Mol failed on {smi!r}: {exc}", file=sys.stderr)
                row = []
            total_smiles += len(row)
            output_smiles.append(row + [""] * (N_OUTPUTS - len(row)))
        log = {
            "start": start_time,
            "end": time.time(),
            "input_smiles": input_smiles,
            "total": total_smiles,
            "expected": N_OUTPUTS * len(input_smiles),
        }
        if is_debug:
            print(
                f"Time taken in seconds: {int(log['end'] - start_time)}, molecules: {total_smiles}"
            )
        return (output_smiles, [], log)
