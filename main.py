import logging

import pandas as pd

from src.models import Experiment, Journal, Ligand, Molecule, Paper
from src.scripts import (create_experiments, create_journals, create_ligands,
                         create_molecules, create_papers, populate_table,
                         upload_mol_props)

pd.options.mode.copy_on_write = True
pd.set_option('future.no_silent_downcasting', True)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    logger.info("Starting...")
    logger.info("Loading dataset")
    df = pd.read_excel("./The_Database_03_12.xlsx")

    logger.info("Creating Journals")
    journals = create_journals()
    logger.info("Uploading Journals to DB")
    populate_table(journals, Journal)

    logger.info("Creating Papers")
    papers = create_papers(df.copy())
    logger.info("Uploading Papers to DB")
    populate_table(papers, Paper)

    logger.info("Creating Molecules")
    molecules = create_molecules(df.copy())
    logger.info("Uploading Molecules to DB")
    populate_table(molecules, Molecule)

    logger.info("Creating Experiments")
    experiments = create_experiments(df.copy())
    logger.info("Uploading Experiments to DB")
    populate_table(experiments, Experiment)

    logger.info("Creating Ligands")
    ligands = create_ligands()
    logger.info("Uploading Ligands")
    populate_table(ligands, Ligand)

    logger.info("Calculate molecules properties")
    upload_mol_props()

    logger.info("Finish.")
