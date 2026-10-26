import os
import random
import sys
import time

import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
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
FP_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def rank_by_similarity(query: str, smiles: "list[str]") -> "list[str]":
    """`smiles` from the most to the least similar to `query`.

    Similarity is the Tanimoto index of Morgan fingerprints (radius 2, 2048 bits); equal
    scores keep their order.
    """
    query_fp = FP_GENERATOR.GetFingerprint(Chem.MolFromSmiles(query))
    fps = [FP_GENERATOR.GetFingerprint(Chem.MolFromSmiles(smile)) for smile in smiles]
    scores = DataStructs.BulkTanimotoSimilarity(query_fp, fps)
    order = sorted(range(len(smiles)), key=scores.__getitem__, reverse=True)
    return [smiles[i] for i in order]


class Mol2MolSimilaritySampler:
    """Samples analogues of a molecule with a similarity-trained Mol2Mol prior of REINVENT 4.

    Each round draws DRAWS_PER_ROUND sequences by multinomial sampling (unseeded, so every call
    gives a different set). Rounds are repeated until N_OUTPUTS distinct molecules that are not the
    input (stereochemistry ignored) are collected, or MAX_ROUNDS is reached. If more than N_OUTPUTS
    were collected, N_OUTPUTS are picked at random. The row is sorted from the most to the least
    similar to the input (Tanimoto) and padded with empty cells.
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
        """Up to N_OUTPUTS distinct analogues of `smiles`, the most similar first."""
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
        picked = random.sample(list(pool), min(N_OUTPUTS, len(pool)))
        return rank_by_similarity(smiles, picked)

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
