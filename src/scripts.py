import re

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Atom, Crippen, Descriptors, Lipinski
from rdkit.Chem.rdchem import MolSanitizeException
from sqlalchemy import Engine, insert, or_, select, text, update

from src.models import Base, Journal, Molecule, Paper
from src.sanitization import sanitize


def populate_table(engine: Engine, df: pd.DataFrame, model: Base, clean_table: bool = False):
    data = df.to_dict(orient="records")

    with engine.connect() as connection:
        if clean_table:
            connection.execute(text(f"DELETE FROM {model.__tablename__}"))

        for rec in data:
            stmt = insert(model).values(**rec)
            connection.execute(stmt)

        connection.commit()


def create_journals() -> pd.DataFrame:
    df = pd.read_excel("Journals_dictionary.xlsx")
    return pd.DataFrame(
        {
            "name": df["New full name"].str.strip(),
            "short_name": df["Shorten name"].str.strip(),
            "abbreviation": df["Abbreviation"].str.strip(),
            "old_name": df["Old name from the database"].str.strip(),
        }
    )


def get_journal_id(journal_name: str, engine: Engine) -> int:
    with engine.connect() as conn:
        stmt = select(Journal.id).where(
            or_(
                Journal.name.ilike(journal_name),
                Journal.short_name.ilike(journal_name),
                Journal.abbreviation.ilike(journal_name),
                Journal.old_name.ilike(journal_name),
            )
        )
        journal_id = conn.scalar(stmt)

        if not journal_id:
            print(f"{journal_name} not found")

    return journal_id


def create_papers(engine: Engine, df: pd.DataFrame) -> pd.DataFrame:
    papers = pd.DataFrame(
        {
            "doi": df["DOI"].str.lower(),
            "title": df["Title"].replace(r"<.*?>", "", regex=True),
            "authors": df["Authors"],
            "abstract": df["Abstract"].str.strip(),
            "journal": df["Journal"].replace(r"&amp;", "&", regex=True),
            "year": df["Year"].astype(int),
            "keywords": df["Keywords"],
            "pages": df["Pages"].replace("Pages Not Available", np.nan),
        }
    )

    papers = papers.drop_duplicates(subset="doi")
    papers["journal_id"] = papers["journal"].apply(get_journal_id, engine=engine).astype(int)
    papers = papers.drop(columns=["journal"])

    return papers


def create_molecules(df: pd.DataFrame) -> pd.DataFrame:
    molecules = pd.DataFrame(
        {
            "cas": df["CAS number"],
            "image_id": df["ImageID"],
            "smiles": df["SMILES"],
            "weight": df["MW (g/mol)"],
        }
    )
    molecules = molecules.drop_duplicates(subset=["cas", "image_id", "smiles"])
    molecules = molecules.round(2)
    molecules["smiles"] = molecules["smiles"].str.replace(
        r"(?<=[A-Za-z])(-|\+)\d*(?!>)", "", regex=True
    )
    return molecules


def get_paper_id(doi: str, engine: Engine) -> int:
    with engine.connect() as conn:
        stmt = select(Paper.id).where(Paper.doi == doi.lower())
        paper_id = conn.scalar(stmt)

        if not paper_id:
            print(f"Did't found paper with doi {doi}")

    return paper_id


def get_molecule_id(row: pd.Series, engine: Engine) -> int:
    cleaned_smiles = re.sub(r"(?<=[A-Za-z])(-|\+)\d*(?!>)", "", row["SMILES"])

    with engine.connect() as conn:
        stmt = (
            select(Molecule.id)
            .where(or_(Molecule.cas == row["CAS number"], Molecule.cas.is_(None)))
            .where(Molecule.image_id == row["ImageID"])
            .where(Molecule.smiles == cleaned_smiles)
        )
        molecule_id = conn.scalar(stmt)

        if not molecule_id:
            print(f"Molecule {row["SMILES"]} not found")

    return molecule_id


def create_source(row: pd.Series) -> pd.Series:
    if pd.notna((row["Homemade"])):
        row["source"] = "homemade"
    elif pd.notna(row["Commercial"]):
        row["source"] = "commercial"
        row["commercial_name"] = row["Commercial"]
    elif pd.notna(row["Clinical"]):
        row["source"] = "clinical"
        row["clinical_name"] = row["Clinical"]
    else:
        row["source"] = pd.NA
        row["commercial_name"] = pd.NA
        row["clinical_name"] = pd.NA
    return row


def freshly_prepared(value: str) -> bool:
    if pd.isna(value):
        return pd.NA
    if value.lower() in ("yes", "freshly prepared"):
        return True
    if value.lower() in ("no", "pre-incubated"):
        return False


def stability_proofs(value: str) -> bool:
    if value is np.nan:
        return pd.NA
    if value == "Y":
        return True
    if value == "N":
        return False


def create_experiments(engine: Engine, df: pd.DataFrame) -> pd.DataFrame:
    experiments = df.drop(
        columns=[
            "ID",
            "Title",
            "Authors",
            "Abstract",
            "Journal",
            "Year",
            "Keywords",
            "Pages",
            "Name",
            "MW (g/mol)",
            "SMILES (UPD)",
        ]
    )
    experiments["paper_id"] = experiments["DOI"].apply(get_paper_id, engine=engine).astype(int)
    experiments = experiments.drop(columns=["DOI"])
    experiments["molecule_id"] = experiments.apply(get_molecule_id, args=(engine,), axis=1).astype(int)
    experiments = experiments.drop(columns=["CAS number", "ImageID", "SMILES"])
    experiments = experiments.replace("no information", np.nan)
    experiments["freshly_prepared"] = experiments["Freshly prepared (yes/no)"].apply(
        freshly_prepared
    )
    experiments = experiments.drop(columns=["Freshly prepared (yes/no)"])
    experiments["stability_proofs"] = experiments["Stability Proofs"].apply(
        stability_proofs
    )
    experiments = experiments.drop(columns=["Stability Proofs"])
    experiments = experiments.apply(create_source, axis=1)
    experiments = experiments.drop(columns=["Homemade", "Commercial", "Clinical"])
    experiments = experiments.melt(
        id_vars=[
            "paper_id",
            "molecule_id",
            "Method",
            "Time",
            "Solvent",
            "freshly_prepared",
            "source",
            "commercial_name",
            "clinical_name",
            "cLogP",
            "LogP",
            "LogD",
            "Solubility (mg/mL)",
            "Log kw",
            "Log k30 (Log k, obtained with a mobile phase containing 30% MeOH)",
            "Log k′0 (Logk, obtained with a mobile phase containing 0% MeOH)",
            "aLogP",
            "LogS",
            "Stability",
            "In vivo",
            "stability_proofs",
            "Stability Proofs Desc",
        ],
        var_name="cell_line",
        value_name="activity",
    )
    experiments = experiments.dropna(subset=["activity"])
    experiments = experiments.drop_duplicates(
        subset=["paper_id", "molecule_id", "Method", "Time", "Solvent", "cell_line"]
    )
    experiments = experiments.rename(
        columns={
            "Method": "method",
            "Time": "time",
            "Solvent": "solvent",
            "LogP": "logP",
            "LogD": "logD",
            "Solubility (mg/mL)": "solubility",
            "Log kw": "log_kw",
            "Log k30 (Log k, obtained with a mobile phase containing 30% MeOH)": "log_k30",
            "Log k′0 (Logk, obtained with a mobile phase containing 0% MeOH)": "log_k0",
            "LogS": "logS",
            "Stability": "stability",
            "In vivo": "in_vivo",
            "Stability Proofs Desc": "stability_proofs_desc",
        }
    )

    return experiments


def generate_ligands(smiles: str, metal_name: str, metal_num: int):
    smiles_wo_isomers = [val for val in smiles.split(".") if metal_name in val]
    ligands_smiles = []
    for val in smiles_wo_isomers:
        try:
            mol = Chem.RWMol(sanitize(val))
        except Exception as e:
            ligands_smiles.append(str(e))
            continue

        pts: list[Atom] = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == metal_num]
        num_of_pts = len(pts)

        for _ in range(num_of_pts):
            pt = pts.pop()

            neigh = pt.GetNeighbors()
            atom: Atom
            for atom in neigh:
                atom.SetBoolProp("metal_bonded", True)

            mol.RemoveAtom(pt.GetIdx())
            pts = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == metal_num]

        ligands = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        for ligand in ligands:
            try:
                Chem.SanitizeMol(ligand)
                ligands_smiles.append((Chem.MolToSmiles(ligand), Chem.MolToCXSmiles(ligand)))
            except MolSanitizeException as e:
                ligands_smiles.append(str(e))
    return ligands_smiles


def create_ligands(engine: Engine, metal_name: str, metal_num: int):
    df = pd.read_sql_table("molecules", engine)
    df = df[["id", "smiles"]]
    df["ligands"] = df["smiles"].apply(generate_ligands, metal_name=metal_name, metal_num=metal_num)
    df = df.drop(columns=["smiles"])
    df = df.explode("ligands", ignore_index=True)
    df = df.dropna(subset=["ligands"])
    df["smiles"] = df["ligands"].apply(lambda x: x[0] if len(x) == 2 else x)
    df["smiles_with_props"] = df["ligands"].apply(lambda x: x[1] if len(x) == 2 else pd.NA)
    df = df.drop(columns=["ligands"])
    df = df.rename(columns={"id": "molecule_id"})
    return df


def calculate_mol_props(row: pd.Series) -> pd.Series:
    try:
        mol = sanitize(row["smiles"])
    except Exception:
        return row

    row["clog_p"] = Crippen.MolLogP(mol)
    row["mr"] = Crippen.MolMR(mol)
    row["nh_oh_count"] = Lipinski.NHOHCount(mol)
    row["no_count"] = Lipinski.NOCount(mol)
    row["num_h_acceptors"] = Lipinski.NumHAcceptors(mol)
    row["num_h_donors"] = Lipinski.NumHDonors(mol)
    row["num_rotatable_bonds"] = Lipinski.NumRotatableBonds(mol)
    row["exact_mol_wt"] = Descriptors.ExactMolWt(mol)
    row["num_of_atoms"] = mol.GetNumAtoms()
    return row


def upload_mol_props(engine: Engine):
    df = pd.read_sql_table("molecules", engine)
    df = df.apply(calculate_mol_props, axis=1)

    data = df.to_dict(orient="records")

    with engine.connect() as connection:
        for rec in data:
            stmt = update(Molecule).where(Molecule.id == rec["id"]).values(**rec)
            connection.execute(stmt)

        connection.commit()
