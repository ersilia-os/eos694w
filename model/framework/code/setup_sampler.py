import os

from torch import device
from reinvent.runmodes import samplers, create_adapter
from reinvent.runmodes.dtos import ChemistryHelpers
from reinvent.chemistry import TransformationTokens


ROOT = os.path.dirname(os.path.abspath(__file__))
CHECKPOINT = os.path.join(ROOT, '..', '..', 'checkpoints')
MODEL = os.path.realpath(os.path.join(CHECKPOINT, 'mol2mol_similarity.prior'))

def setup_sampler(batch_size: int, chemistry: ChemistryHelpers) -> samplers.Mol2MolSampler:
  """Setup the sampling module.
  
  Args:
  -----

  batch_size: int
    Number of molecules to generate.

  """
  temperature = 1.0
  unique_sequences = True
  isomeric = True
  sample_strategy = "beamsearch"
  tokens = TransformationTokens()
  randomize_smiles = False # For Mol2Mol reinvent 4 sets randomize_smiles to False.


  agent, _, _ = create_adapter(MODEL, "inference", device('cpu'))

  sampler =  samplers.Mol2MolSampler(
    agent,
    batch_size=batch_size,
    sample_strategy=sample_strategy,
    isomeric=isomeric,
    randomize_smiles=randomize_smiles,
    unique_sequences=unique_sequences,
    chemistry=chemistry,
    temperature=temperature,
    tokens=tokens
  )

  return sampler

