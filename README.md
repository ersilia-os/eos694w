# REINVENT 4 Mol2MolMediumSimilarity

The Mol2Mol Medium Similarity option of REINVENT4 generates a diverse set of approximately 100 novel small molecules. These generated molecules are designed to retain a moderate level of structural and chemical similarity to the input compound, enabling scaffold exploration while preserving key pharmacophoric or functional features. This approach is particularly useful for hit expansion, analog generation, and chemical space exploration tasks in early-stage drug discovery workflows.

This model was incorporated on 2024-02-07.Last packaged on 2026-02-05.

## Information
### Identifiers
- **Ersilia Identifier:** `eos694w`
- **Slug:** `reinvent4-mol2mol-medium-similarity`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Generation`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Similarity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `100`
- **Output Consistency:** `Variable`
- **Interpretation:** Model generates up to 100 similar molecules per input molecule.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smi_00 | string |  | Generated molecule index 0 using the mol2mol medium similarity prior from REINVENT |
| smi_01 | string |  | Generated molecule index 1 using the mol2mol medium similarity prior from REINVENT |
| smi_02 | string |  | Generated molecule index 2 using the mol2mol medium similarity prior from REINVENT |
| smi_03 | string |  | Generated molecule index 3 using the mol2mol medium similarity prior from REINVENT |
| smi_04 | string |  | Generated molecule index 4 using the mol2mol medium similarity prior from REINVENT |
| smi_05 | string |  | Generated molecule index 5 using the mol2mol medium similarity prior from REINVENT |
| smi_06 | string |  | Generated molecule index 6 using the mol2mol medium similarity prior from REINVENT |
| smi_07 | string |  | Generated molecule index 7 using the mol2mol medium similarity prior from REINVENT |
| smi_08 | string |  | Generated molecule index 8 using the mol2mol medium similarity prior from REINVENT |
| smi_09 | string |  | Generated molecule index 9 using the mol2mol medium similarity prior from REINVENT |

_10 of 100 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos694w](https://hub.docker.com/r/ersiliaos/eos694w)
- **Docker Architecture:** `AMD64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos694w.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos694w.zip)

### Resource Consumption
- **Model Size (Mb):** `233`
- **Environment Size (Mb):** `2360`
- **Image Size (Mb):** `2554.15`

**Computational Performance (seconds):**
- 10 inputs: `86.12`
- 100 inputs: `-1`
- 10000 inputs: `-1`

### References
- **Source Code**: [https://github.com/MolecularAI/REINVENT4](https://github.com/MolecularAI/REINVENT4)
- **Publication**: [https://chemrxiv.org/engage/chemrxiv/article-details/65463cafc573f893f1cae33a](https://chemrxiv.org/engage/chemrxiv/article-details/65463cafc573f893f1cae33a)
- **Publication Type:** `Preprint`
- **Publication Year:** `2023`
- **Ersilia Contributor:** [ankitskvmdam](https://github.com/ankitskvmdam)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [Apache-2.0](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos694w
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos694w
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
