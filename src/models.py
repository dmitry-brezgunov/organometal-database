from decimal import Decimal
from typing import Optional

from sqlalchemy import ForeignKey, Numeric, UniqueConstraint
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    id: Mapped[int] = mapped_column(primary_key=True, sort_order=-1)


class Journal(Base):
    __tablename__ = "journals"
    name: Mapped[str] = mapped_column(unique=True)
    short_name: Mapped[str] = mapped_column(unique=True)
    abbreviation: Mapped[str]
    old_name: Mapped[str] = mapped_column(unique=True)
    papers: Mapped[set["Paper"]] = relationship(back_populates="journal")

    def __repr__(self) -> str:
        return f"Journal(id={self.id!r}, name={self.name!r})"


class Paper(Base):
    __tablename__ = "papers"
    doi: Mapped[str] = mapped_column(unique=True)
    title: Mapped[str] = mapped_column(unique=True)
    authors: Mapped[str]
    abstract: Mapped[str]
    journal_id: Mapped[int] = mapped_column(ForeignKey("journals.id"))
    journal: Mapped["Journal"] = relationship(back_populates="papers")
    year: Mapped[int]
    keywords: Mapped[Optional[str]]
    pages: Mapped[Optional[str]]
    experiments: Mapped[set["Experiment"]] = relationship(back_populates="paper")


class Molecule(Base):
    __tablename__ = "molecules"
    __table_args__ = (
        UniqueConstraint(
            "cas", "image_id", "smiles", name="cas_image_smiles_constraint"
        ),
    )
    cas: Mapped[Optional[str]]
    image_id: Mapped[str]
    smiles: Mapped[str]
    weight: Mapped[Optional[Decimal]] = mapped_column(Numeric(6, 2))
    experiments: Mapped[set["Experiment"]] = relationship(back_populates="molecule")
    ligands: Mapped[set["Ligand"]] = relationship(back_populates="molecule")
    clog_p: Mapped[Optional[float]]
    mr: Mapped[Optional[float]]
    nh_oh_count: Mapped[Optional[int]]
    no_count: Mapped[Optional[int]]
    num_h_acceptors: Mapped[Optional[int]]
    num_h_donors: Mapped[Optional[int]]
    num_rotatable_bonds: Mapped[Optional[int]]
    exact_mol_wt: Mapped[Optional[float]]
    num_of_atoms: Mapped[Optional[int]]
    name: Mapped[Optional[str]]


class Experiment(Base):
    __tablename__ = "experiments"
    __table_args__ = (
        UniqueConstraint(
            "paper_id",
            "molecule_id",
            "method",
            "time",
            "solvent",
            "cell_line",
            name="experiments_constraint",
        ),
    )
    paper_id: Mapped[int] = mapped_column(ForeignKey("papers.id"))
    paper: Mapped["Paper"] = relationship(back_populates="experiments")
    molecule_id: Mapped[int] = mapped_column(ForeignKey("molecules.id"))
    molecule: Mapped["Molecule"] = relationship(back_populates="experiments")
    method: Mapped[Optional[str]]
    time: Mapped[Optional[Decimal]] = mapped_column(Numeric(4, 1))
    solvent: Mapped[Optional[str]]
    freshly_prepared: Mapped[Optional[bool]]
    source: Mapped[Optional[str]]
    commercial_name: Mapped[Optional[str]]
    clinical_name: Mapped[Optional[str]]
    cLogP: Mapped[Optional[str]]
    logP: Mapped[Optional[str]]
    logD: Mapped[Optional[str]]
    solubility: Mapped[Optional[str]]
    log_kw: Mapped[Optional[str]]
    log_k30: Mapped[Optional[str]]
    log_k0: Mapped[Optional[str]]
    aLogP: Mapped[Optional[str]]
    logS: Mapped[Optional[str]]
    stability: Mapped[Optional[str]]
    in_vivo: Mapped[Optional[str]]
    stability_proofs: Mapped[Optional[bool]]
    stability_proofs_desc: Mapped[Optional[str]]
    cell_line: Mapped[str]
    activity: Mapped[str]


class Ligand(Base):
    __tablename__ = "ligands"
    molecule_id: Mapped[int] = mapped_column(ForeignKey("molecules.id"))
    molecule: Mapped["Molecule"] = relationship(back_populates="ligands")
    smiles: Mapped[str]
    smiles_with_props: Mapped[Optional[str]]
