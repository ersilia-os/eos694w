# imports
import os
import sys

import torch
from reinvent.runmodes.samplers.run_sampling import run_sampling


# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# my model
def my_model():
     
    # Path of the model.
    MODEL = os.path.realpath(
      os.path.join(
        os.path.dirname(__file__),
        '..',
        '..',
        'checkpoints',
        'prior',
        'mol2mol_similarity.prior'
        )
      )

    # Default config
    config = {
        "run_type": "sampling",
        "parameters": {
            "model_file": MODEL,
            "unique_molecules": True,
            "num_smiles": 100,
            "output_file": output_file,
            "smiles_file": input_file,
            "sample_strategy": "beamsearch",
            "temperature": 1.0
        }
    }

    run_sampling(config, torch.device('cpu'))



# run model
my_model()