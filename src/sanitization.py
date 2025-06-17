import re

from rdkit import Chem

DATIVE_VALENCE = {
    6: [2, 4],
    7: [3],
    8: [2],
    15: [3, 5],
    16: [2, 4, 6],
    33: [3, 5],
    34: [2, 4, 6]
}

NBP_CHARGES = {5: (4, -1), 7: (4, 1), 15: (6, -1)}


def is_transition_metal(at: Chem.Atom) -> bool:
    n = at.GetAtomicNum()
    return (n >= 21 and n <= 31) or (n >= 39 and n <= 48) or (n >= 71 and n <= 80)


def is_metal(at: Chem.Atom) -> bool:
    metals = (3, 4, 11, 12, 13, 19, 20, 31, 37, 38, 49, 50, 55, 56, 81, 82, 83, 87, 88)
    n = at.GetAtomicNum()
    return (n in metals) or is_transition_metal(at)


def get_max_valence(at_num: int, actual_val: int) -> int:
    valence_list = DATIVE_VALENCE[at_num]

    if len(valence_list) == 1:
        return valence_list[0]
    if actual_val in valence_list:
        return actual_val
    if actual_val - 1 in valence_list:
        return actual_val - 1
    if actual_val <= valence_list[0]:
        return valence_list[0]
    if actual_val >= valence_list[-1]:
        return valence_list[-1]


def set_dative_bonds(mol: Chem.Mol) -> Chem.RWMol:
    """convert some bonds to dative

    Replaces some single bonds between metals and atoms with atomic numbers in fomAtoms
    with dative bonds. The replacement is only done if the atom has "too many" bonds.

    Returns the modified molecule.

    """
    rwmol = Chem.RWMol(mol)
    rwmol.UpdatePropertyCache(strict=False)
    metals = [at for at in rwmol.GetAtoms() if is_transition_metal(at)]
    for metal in metals:
        dative_nbr = [at for at in metal.GetNeighbors() if at.GetAtomicNum() in DATIVE_VALENCE]
        for nbr in dative_nbr:
            nbr_atom_num = nbr.GetAtomicNum()
            nbr_actual_valence = nbr.GetExplicitValence()
            nbr_max_valence = get_max_valence(nbr_atom_num, nbr_actual_valence)
            bond_type = rwmol.GetBondBetweenAtoms(nbr.GetIdx(), metal.GetIdx()).GetBondType()
            if nbr_actual_valence > nbr_max_valence and bond_type == Chem.BondType.SINGLE:
                rwmol.RemoveBond(nbr.GetIdx(), metal.GetIdx())
                rwmol.AddBond(nbr.GetIdx(), metal.GetIdx(), Chem.BondType.DATIVE)
    return rwmol


def set_charges_for_nbp(mol: Chem.Mol) -> Chem.Mol:
    mol.UpdatePropertyCache(strict=False)
    atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in (5, 7, 15)]
    for atom in atoms:
        metals_nbr = [at for at in atom.GetNeighbors() if is_metal(at)]
        num_bonds = atom.GetExplicitValence()
        atom_num = atom.GetAtomicNum()

        if not metals_nbr and num_bonds == NBP_CHARGES[atom_num][0]:
            atom.SetFormalCharge(NBP_CHARGES[atom_num][1])

    return mol


def set_charges_for_pt(mol: Chem.Mol) -> Chem.Mol:
    mol.UpdatePropertyCache(strict=False)
    pts = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 78]

    for atom in pts:
        num_bonds = len(atom.GetBonds())
        num_non_dative_bonds = len([bond for bond in atom.GetBonds() if bond.GetBondType() != Chem.BondType.DATIVE])

        if num_bonds == 6:
            atom.SetFormalCharge(4 - num_non_dative_bonds)
        elif num_bonds == 4:
            atom.SetFormalCharge(2 - num_non_dative_bonds)

    return mol


def sanitize(mol: str | Chem.Mol | Chem.RWMol) -> Chem.Mol | Chem.RWMol:
    if isinstance(mol, str):
        cleaned_smiles = re.sub(r"(?<=[A-Za-z])(-|\+)\d*(?!>)", "", mol)
        mol = Chem.MolFromSmiles(cleaned_smiles, sanitize=False)

    mol = set_dative_bonds(mol)
    mol = set_charges_for_nbp(mol)
    Chem.SanitizeMol(mol)
    mol = set_charges_for_pt(mol)
    return mol
