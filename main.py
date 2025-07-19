import logging

import pandas as pd
from sqlalchemy import create_engine

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
    df = pd.read_excel("./cuprum/full_db.xlsx")
    engine = create_engine("sqlite:///cuprum.db")

    logger.info("Creating Journals")
    journals = create_journals()
    logger.info("Uploading Journals to DB")
    populate_table(engine, journals, Journal, clean_table=True)

    logger.info("Creating Papers")
    papers = create_papers(engine, df.copy())
    logger.info("Uploading Papers to DB")
    populate_table(engine, papers, Paper, clean_table=True)

    logger.info("Creating Molecules")
    molecules = create_molecules(df.copy())
    logger.info("Uploading Molecules to DB")
    populate_table(engine, molecules, Molecule, clean_table=True)

    logger.info("Creating Experiments")
    experiments = create_experiments(engine, df.copy())
    logger.info("Uploading Experiments to DB")
    populate_table(engine, experiments, Experiment, clean_table=True)

    logger.info("Creating Ligands")
    ligands = create_ligands(engine, metal_name="Cu", metal_num=29)
    logger.info("Uploading Ligands")
    populate_table(engine, ligands, Ligand, clean_table=True)

    logger.info("Calculate molecules properties")
    upload_mol_props(engine)

    logger.info("Finish.")
