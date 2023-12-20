#!/usr/bin/python
# -*- coding: latin-1 -*-
from Bio import PDB
from Bio.PDB import PDBIO
from Bio.SVDSuperimposer import SVDSuperimposer
import Bio.PDB
import numpy 


def calculate_rmsd(atoms1, atoms2):
    """
    Oblicza miarê RMSD dla dwóch zestawów atomów.
    """
    assert len(atoms1) == len(atoms2)

    N = len(atoms1)
    sum_sq_dist = 0.0

    for i in range(N):
        dist = numpy.linalg.norm(atoms1[i] - atoms2[i])
        sum_sq_dist += dist ** 2

    rmsd = numpy.sqrt(sum_sq_dist / N)
    return rmsd


pdb_parser = PDB.PDBParser(QUIET=True)
structure = pdb_parser.get_structure("R1107_reference", "C:/Users/anna/Desktop/pp/bioinf_strukturalna/zad5/R1107/R1107_reference.pdb")
ref_model = structure[0]

model_filename = "C:/Users/anna/Desktop/pp/bioinf_strukturalna/zad5/R1107/R1107TS081.pdb"
models_structure = PDB.PDBParser(QUIET=True).get_structure("models", model_filename)
num_models = len(models_structure)

rmsd_results = []

for i, alt_model in enumerate(models_structure):
    ref_atoms = []
    alt_atoms = []

    for (ref_chain, alt_chain) in zip(ref_model, alt_model):
        for ref_res, alt_res in zip(ref_chain, alt_chain):
            for ref_atom, alt_atom in zip(ref_res, alt_res):
                if ref_atom.name == alt_atom.name:
                    ref_atoms.append(ref_atom)
                    alt_atoms.append(alt_atom)


    superimposer = PDB.Superimposer()
    superimposer.set_atoms(ref_atoms, alt_atoms)

    superimposer.apply(alt_model.get_atoms())

    rmsd = calculate_rmsd(ref_atoms, alt_atoms)
    rmsd_results.append((i+1, rmsd))

    print(f"RMSD dla Modelu {i+1}: {rmsd} {superimposer.rms}")


best_model, best_rmsd = min(rmsd_results, key=lambda x: x[1])
worst_model, worst_rmsd = max(rmsd_results, key=lambda x: x[1])


best_model_structure = models_structure[best_model-1]
worst_model_structure = models_structure[worst_model-1]

io = PDBIO()
io.set_structure(best_model_structure)
io.save(f"best_model.pdb", write_end=True)

io.set_structure(worst_model_structure)
io.save(f"worst_model.pdb", write_end=True)

print(f"Najlepszy model: {best_model}, RMSD: {best_rmsd}")
print(f"Najgorszy model: {worst_model}, RMSD: {worst_rmsd}")





