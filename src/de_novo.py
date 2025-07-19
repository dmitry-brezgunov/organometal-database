import itertools
import warnings
from math import factorial

import pandas as pd
from rdkit import Chem, RDLogger
from sqlalchemy import Engine
from tqdm import tqdm

from src.sanitization import sanitize

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")


def generate_de_novo_from_ligands(engine: Engine, query: str):
    df = pd.read_sql(query, engine)
    ligands = df["smiles_with_props"].to_list()

    smiles_molecules = dict()
    repeats = 6
    length = factorial(len(ligands) + repeats - 1) / factorial(repeats) / factorial(len(ligands) - 1)
    for combination in tqdm(itertools.combinations_with_replacement(ligands, 6), total=int(length)):
        mol = Chem.Mol()
        coord_num = 0
        for ligand_smiles in combination:
            ligand = Chem.MolFromSmiles(ligand_smiles)
            dentat = len([atom for atom in ligand.GetAtoms() if atom.HasProp("metal_bonded")])
            coord_num += dentat
            mol = Chem.CombineMols(mol, ligand)

            if coord_num in (4, 6) and Chem.MolToSmiles(mol) not in smiles_molecules:
                smiles_molecules[Chem.MolToSmiles(mol)] = mol

            if coord_num >= 6:
                break

    de_novo = set()
    for mol in smiles_molecules.values():
        rwmol = Chem.RWMol(mol)
        pt_idx = rwmol.AddAtom(Chem.Atom(78))
        for atom in rwmol.GetAtoms():
            if atom.HasProp("metal_bonded"):
                rwmol.AddBond(atom.GetIdx(), pt_idx, Chem.BondType.SINGLE)

        try:
            rwmol = sanitize(rwmol)
            de_novo.add(Chem.MolToSmiles(rwmol))
        except Chem.MolSanitizeException:
            continue

    organometaldb = pd.read_sql("SELECT DISTINCT smiles FROM molecules", engine)
    scifinder = pd.read_csv("scifinder-csd-wo-dative.csv")
    existing_smiles = pd.concat([organometaldb, scifinder])
    existing_smiles = existing_smiles.drop_duplicates(subset=["smiles"])

    for smiles in tqdm(existing_smiles["smiles"].to_list()):
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        try:
            mol = sanitize(mol)
            canonical_smiles = Chem.MolToSmiles(mol)
            de_novo.discard(canonical_smiles)
        except (Chem.MolSanitizeException, ValueError):
            continue

    print(f"Total denovo: {len(smiles_molecules)}")
    print(f"Valid unique de novo: {len(de_novo)}")

    pd.DataFrame({"smiles": list(de_novo)}).to_csv("de_novo_smart.csv", index=False)


if __name__ == "__main__":
    query = """
        WITH cleaned_activity as (
            SELECT
                m.id,
                m.smiles,
                CASE
                    WHEN INSTR(e.activity, '±') THEN CAST(SUBSTR(e.activity, 1, INSTR(e.activity, '±') - 1) as float)
                    WHEN INSTR(e.activity, '<') THEN CAST(SUBSTR(e.activity, INSTR(e.activity, '<') + 1) as float)
                    WHEN INSTR(e.activity, '>') THEN CAST(SUBSTR(e.activity, INSTR(e.activity, '>') + 1) as float)
                    WHEN e.activity = 'n.d.' THEN 9999
                    ELSE CAST(e.activity as float)
                END as activity,
                clog_p,
                mr,
                nh_oh_count,
                no_count,
                num_h_acceptors,
                num_h_donors,
                num_rotatable_bonds,
                exact_mol_wt,
                num_of_atoms
            FROM experiments e
            JOIN molecules m ON e.molecule_id = m.id
            WHERE e.cell_line in ('MCF-7', 'A549', 'A2780', 'A2780cis') and e.time = 72
        ),
        molecules_filtered as (
            select
                id,
                smiles,
                max(activity) as activity,
                clog_p,
                mr,
                nh_oh_count,
                no_count,
                num_h_acceptors,
                num_h_donors,
                num_rotatable_bonds,
                exact_mol_wt,
                num_of_atoms
            from cleaned_activity
            where
                activity < 5
                and num_h_donors <= 5
                and num_h_acceptors <= 10
                and exact_mol_wt < 500
                and clog_p <= 5
                and mr between 40 and 130
                and num_of_atoms between 20 and 70
                and num_rotatable_bonds <= 10
            group by smiles
            order by activity
        )
        select smiles, smiles_with_props, count(molecule_id) as cnt
        from ligands
        where molecule_id in (select id from molecules_filtered)
        group by smiles, smiles_with_props
        order by cnt desc
        limit 20
    """
    generate_de_novo_from_ligands(query)
