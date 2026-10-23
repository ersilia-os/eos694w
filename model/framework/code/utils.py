from rdkit import Chem


def is_same_2d_structure(smiles1: str, smiles2: str) -> bool:
    """True if both SMILES are the same molecule once stereochemistry is ignored.

    The canonical SMILES without stereochemistry are compared. REINVENT 4 writes the
    input in a different form and returns stereoisomers of it, so comparing strings,
    or fingerprints, does not catch every copy of the input.
    """
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return False
    return Chem.MolToSmiles(mol1, isomericSmiles=False) == Chem.MolToSmiles(
        mol2, isomericSmiles=False
    )
