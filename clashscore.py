#!/usr/bin/python
# -*- coding: latin-1 -*-
from Bio import PDB
import numpy as np


def epsilon_vdw_PDB(atom_name, residue_name) :

    dR = {}
    dR["C"] = 1.9080
    dR["CA"] = 1.9080
    dR["CM"] = 1.9080
    dR["Cs"] = 3.3950
    dR["CT"] = 1.9080
    dR["H"] = 0.6000
    dR["H1"] = 1.3870
    dR["H2"] = 1.2870
    dR["H3"] = 1.1870
    dR["H4"] = 1.4090
    dR["H5"] = 1.3590
    dR["HA"] = 1.4590
    dR["HC"] = 1.4870
    dR["HO"] = 0.0000
    dR["HP"] = 1.1000
    dR["HS"] = 0.6000
    dR["HW"] = 0.0000
    dR["IP"] = 1.8680
    dR["K"] = 2.6580
    dR["Li"] = 1.1370
    dR["N2m"] = 1.8240
    dR["N3n"] = 1.8750
    dR["O"] = 1.6612
    dR["O2"] = 1.6612
    dR["OH"] = 1.7210
    dR["OS"] = 1.6837
    dR["OW"] = 1.7683
    dR["S"] = 2.0000
    dR["SH"] = 2.0000

    dE = {}
    dE["C"] = 0.0860
    dE["CA"] = 0.0860
    dE["CM"] = 0.0860
    dE["Cs"] = 0.0000806
    dE["CT"] = 0.1094
    dE["H"] = 0.0157
    dE["H1"] = 0.0157
    dE["H2"] = 0.0157
    dE["H3"] = 0.0157
    dE["H4"] = 0.0150
    dE["H5"] = 0.0150
    dE["HA"] = 0.0150
    dE["HC"] = 0.0157
    dE["HO"] = 0.0000
    dE["HP"] = 0.0157
    dE["HS"] = 0.0157
    dE["HW"] = 0.0000
    dE["IP"] = 0.00277
    dE["K"] = 0.000328
    dE["Li"] = 0.0183
    dE["N2m"] = 0.1700
    dE["N3n"] = 0.1700
    dE["O"] = 0.2100
    dE["O2"] = 0.2100
    dE["OH"] = 0.2104
    dE["OS"] = 0.1700
    dE["OW"] = 0.1520
    dE["S"] = 0.2500
    dE["SH"] = 0.2500

    # vdw radius

    dvdw = {}
    dvdw["GLY"] = {}
    dvdw["GLY"]["N"] = dR["N2m"]
    dvdw["GLY"]["CA"] = dR["CT"]
    dvdw["GLY"]["C"] = dR["C"]
    dvdw["GLY"]["O"] = dR["O"]
    dvdw["GLY"]["H"] = dR["H"]
    dvdw["GLY"]["HA2"] = dR["H1"]
    dvdw["GLY"]["HA3"] = dR["H1"]
    
    dvdw["ALA"] = {}
    dvdw["ALA"]["N"] = dR["N2m"]
    dvdw["ALA"]["CA"] = dR["CT"]
    dvdw["ALA"]["C"] = dR["C"]
    dvdw["ALA"]["O"] = dR["O"]
    dvdw["ALA"]["CB"] = dR["CT"]
    dvdw["ALA"]["H"] = dR["H"]
    dvdw["ALA"]["HB1"] = dR["HC"]
    dvdw["ALA"]["HB2"] = dR["HC"]
    dvdw["ALA"]["HB3"] = dR["HC"]
    dvdw["ALA"]["HA"] = dR["H1"]

    dvdw["ASP"] = {}
    dvdw["ASP"]["N"] = dR["N2m"]
    dvdw["ASP"]["CA"] = dR["CT"]
    dvdw["ASP"]["C"] = dR["C"]
    dvdw["ASP"]["O"] = dR["O"]
    dvdw["ASP"]["CB"] = dR["CT"]
    dvdw["ASP"]["CG"] = dR["C"]
    dvdw["ASP"]["OD1"] = dR["O2"]
    dvdw["ASP"]["OD2"] = dR["O2"]
    dvdw["ASP"]["H"] = dR["H"]
    dvdw["ASP"]["HA"] = dR["H1"]
    dvdw["ASP"]["HB2"] = dR["HC"]
    dvdw["ASP"]["HB3"] = dR["HC"]

    dvdw["GLU"] = {}
    dvdw["GLU"]["N"] = dR["N2m"]
    dvdw["GLU"]["CA"] = dR["CT"]
    dvdw["GLU"]["C"] = dR["C"]
    dvdw["GLU"]["O"] = dR["O"]
    dvdw["GLU"]["CB"] = dR["CT"]
    dvdw["GLU"]["CG"] = dR["CT"]
    dvdw["GLU"]["CD"] = dR["C"]
    dvdw["GLU"]["OE1"] = dR["O2"]
    dvdw["GLU"]["OE2"] = dR["O2"]
    dvdw["GLU"]["H"] = dR["H"]
    dvdw["GLU"]["HA"] = dR["H1"]
    dvdw["GLU"]["HB2"] = dR["HC"]
    dvdw["GLU"]["HB3"] = dR["HC"]
    dvdw["GLU"]["HG2"] = dR["HC"]
    dvdw["GLU"]["HG3"] = dR["HC"]

    dvdw["LEU"] = {}
    dvdw["LEU"]["N"] = dR["N2m"]
    dvdw["LEU"]["CA"] = dR["CT"]
    dvdw["LEU"]["C"] = dR["C"]
    dvdw["LEU"]["O"] = dR["O"]
    dvdw["LEU"]["CB"] = dR["CT"]
    dvdw["LEU"]["CG"] = dR["CT"]
    dvdw["LEU"]["CD1"] = dR["CT"]
    dvdw["LEU"]["CD2"] = dR["CT"]
    dvdw["LEU"]["H"] = dR["H"]
    dvdw["LEU"]["HA"] = dR["H1"]
    dvdw["LEU"]["HB2"] = dR["HC"]
    dvdw["LEU"]["HB3"] = dR["HC"]
    dvdw["LEU"]["HG"] = dR["HC"]
    dvdw["LEU"]["HD11"] = dR["HC"]
    dvdw["LEU"]["HD12"] = dR["HC"]
    dvdw["LEU"]["HD13"] = dR["HC"]
    dvdw["LEU"]["HD21"] = dR["HC"]
    dvdw["LEU"]["HD22"] = dR["HC"]
    dvdw["LEU"]["HD23"] = dR["HC"]
    
    dvdw["ASN"] = {}
    dvdw["ASN"]["N"] = dR["N2m"]
    dvdw["ASN"]["CA"] = dR["CT"]
    dvdw["ASN"]["C"] = dR["C"]
    dvdw["ASN"]["O"] = dR["O"]
    dvdw["ASN"]["CB"] = dR["CT"]
    dvdw["ASN"]["CG"] = dR["C"]
    dvdw["ASN"]["ND2"] = dR["N2m"]
    dvdw["ASN"]["OD1"] = dR["O"]
    dvdw["ASN"]["H"] = dR["H"]
    dvdw["ASN"]["HA"] = dR["H1"]
    dvdw["ASN"]["HB2"] = dR["HC"]
    dvdw["ASN"]["HB3"] = dR["HC"]
    dvdw["ASN"]["HD21"] = dR["H"]
    dvdw["ASN"]["HD22"] = dR["H"]

    dvdw["GLN"] = {}
    dvdw["GLN"]["N"] = dR["N2m"]
    dvdw["GLN"]["CA"] = dR["CT"]
    dvdw["GLN"]["C"] = dR["C"]
    dvdw["GLN"]["O"] = dR["O"]
    dvdw["GLN"]["CB"] = dR["CT"]
    dvdw["GLN"]["CG"] = dR["CT"]
    dvdw["GLN"]["CD"] = dR["C"]
    dvdw["GLN"]["NE2"] = dR["N2m"]
    dvdw["GLN"]["OE1"] = dR["O"]
    dvdw["GLN"]["H"] = dR["H"]
    dvdw["GLN"]["HA"] = dR["H1"]
    dvdw["GLN"]["HB2"] = dR["HC"]
    dvdw["GLN"]["HB3"] = dR["HC"]
    dvdw["GLN"]["HG2"] = dR["HC"]
    dvdw["GLN"]["HG3"] = dR["HC"]
    dvdw["GLN"]["HE21"] = dR["H"]
    dvdw["GLN"]["HE22"] = dR["H"]

    dvdw["ILE"] = {}
    dvdw["ILE"]["N"] = dR["N2m"]
    dvdw["ILE"]["CA"] = dR["CT"]
    dvdw["ILE"]["C"] = dR["C"]
    dvdw["ILE"]["O"] = dR["O"]
    dvdw["ILE"]["CB"] = dR["CT"]
    dvdw["ILE"]["CG1"] = dR["CT"]
    dvdw["ILE"]["CG2"] = dR["CT"]
    dvdw["ILE"]["CD1"] = dR["CT"]
    dvdw["ILE"]["H"] = dR["H"]
    dvdw["ILE"]["HA"] = dR["H1"]
    dvdw["ILE"]["HB"] = dR["HC"]
    dvdw["ILE"]["HG12"] = dR["HC"]
    dvdw["ILE"]["HG13"] = dR["HC"]
    dvdw["ILE"]["HG21"] = dR["HC"]
    dvdw["ILE"]["HG22"] = dR["HC"]
    dvdw["ILE"]["HG23"] = dR["HC"]
    dvdw["ILE"]["HD11"] = dR["HC"]
    dvdw["ILE"]["HD12"] = dR["HC"]
    dvdw["ILE"]["HD13"] = dR["HC"]

    dvdw["VAL"] = {}
    dvdw["VAL"]["N"] = dR["N2m"]
    dvdw["VAL"]["CA"] = dR["CT"]
    dvdw["VAL"]["C"] = dR["C"]
    dvdw["VAL"]["O"] = dR["O"]
    dvdw["VAL"]["CB"] = dR["CT"]
    dvdw["VAL"]["CG1"] = dR["CT"]
    dvdw["VAL"]["CG2"] = dR["CT"]
    dvdw["VAL"]["H"] = dR["H"]
    dvdw["VAL"]["HA"] = dR["H1"]
    dvdw["VAL"]["HB"] = dR["HC"]
    dvdw["VAL"]["HG11"] = dR["HC"]
    dvdw["VAL"]["HG12"] = dR["HC"]
    dvdw["VAL"]["HG13"] = dR["HC"]
    dvdw["VAL"]["HG21"] = dR["HC"]
    dvdw["VAL"]["HG22"] = dR["HC"]
    dvdw["VAL"]["HG23"] = dR["HC"]

    dvdw["SER"] = {}
    dvdw["SER"]["N"] = dR["N2m"]
    dvdw["SER"]["CA"] = dR["CT"]
    dvdw["SER"]["C"] = dR["C"]
    dvdw["SER"]["O"] = dR["O"]
    dvdw["SER"]["CB"] = dR["CT"]
    dvdw["SER"]["OG"] = dR["OH"]
    dvdw["SER"]["H"] = dR["H"]
    dvdw["SER"]["HA"] = dR["H1"]
    dvdw["SER"]["HB2"] = dR["H1"]
    dvdw["SER"]["HB3"] = dR["H1"]
    dvdw["SER"]["HG"] = dR["HO"]
 
    dvdw["THR"] = {}
    dvdw["THR"]["N"] = dR["N2m"]
    dvdw["THR"]["CA"] = dR["CT"]
    dvdw["THR"]["C"] = dR["C"]
    dvdw["THR"]["O"] = dR["O"]
    dvdw["THR"]["CB"] = dR["CT"]
    dvdw["THR"]["OG1"] = dR["OH"]
    dvdw["THR"]["CG2"] = dR["CT"]
    dvdw["THR"]["H"] = dR["H"]
    dvdw["THR"]["HA"] = dR["H1"]
    dvdw["THR"]["HB"] = dR["H1"]
    dvdw["THR"]["HG1"] = dR["HO"]
    dvdw["THR"]["HG21"] = dR["HC"]
    dvdw["THR"]["HG22"] = dR["HC"]
    dvdw["THR"]["HG23"] = dR["HC"]

    dvdw["CYS"] = {}
    dvdw["CYS"]["N"] = dR["N2m"]
    dvdw["CYS"]["CA"] = dR["CT"]
    dvdw["CYS"]["C"] = dR["C"]
    dvdw["CYS"]["O"] = dR["O"]
    dvdw["CYS"]["CB"] = dR["CT"]
    dvdw["CYS"]["SG"] = dR["SH"]
    dvdw["CYS"]["H"] = dR["H"]
    dvdw["CYS"]["HA"] = dR["H1"]
    dvdw["CYS"]["HB1"] = dR["H1"]
    dvdw["CYS"]["HB2"] = dR["H1"]
    dvdw["CYS"]["HG"] = dR["HS"]

    dvdw["PRO"] = {}
    dvdw["PRO"]["N"] = dR["N2m"]
    dvdw["PRO"]["CA"] = dR["CT"]
    dvdw["PRO"]["C"] = dR["C"]
    dvdw["PRO"]["O"] = dR["O"]
    dvdw["PRO"]["CB"] = dR["CT"]
    dvdw["PRO"]["CG"] = dR["CT"]
    dvdw["PRO"]["CD"] = dR["CT"]
    dvdw["PRO"]["HA"] = dR["H1"]
    dvdw["PRO"]["HB2"] = dR["HC"]
    dvdw["PRO"]["HB3"] = dR["HC"]
    dvdw["PRO"]["HG2"] = dR["HC"]
    dvdw["PRO"]["HG3"] = dR["HC"]
    dvdw["PRO"]["HD2"] = dR["H1"]
    dvdw["PRO"]["HD3"] = dR["H1"]
    
    dvdw["ARG"] = {}
    dvdw["ARG"]["N"] = dR["N2m"]
    dvdw["ARG"]["CA"] = dR["CT"]
    dvdw["ARG"]["C"] = dR["C"]
    dvdw["ARG"]["O"] = dR["O"]
    dvdw["ARG"]["CB"] = dR["CT"]
    dvdw["ARG"]["CG"] = dR["CT"]
    dvdw["ARG"]["CD"] = dR["CT"]
    dvdw["ARG"]["NE"] = dR["N2m"]
    dvdw["ARG"]["CZ"] = dR["CA"]
    dvdw["ARG"]["NH1"] = dR["N2m"]
    dvdw["ARG"]["NH2"] = dR["N2m"]
    dvdw["ARG"]["H"] = dR["H"]
    dvdw["ARG"]["HA"] = dR["H1"]
    dvdw["ARG"]["HB2"] = dR["HC"]
    dvdw["ARG"]["HB3"] = dR["HC"]
    dvdw["ARG"]["HG2"] = dR["HC"]
    dvdw["ARG"]["HG3"] = dR["HC"]
    dvdw["ARG"]["HD2"] = dR["H1"]
    dvdw["ARG"]["HD3"] = dR["H1"]
    dvdw["ARG"]["HE"] = dR["H"]
    dvdw["ARG"]["HH11"] = dR["H"]
    dvdw["ARG"]["HH12"] = dR["H"]
    dvdw["ARG"]["HH21"] = dR["H"]  
    dvdw["ARG"]["HH22"] = dR["H"]

    dvdw["LYS"] = {}
    dvdw["LYS"]["N"] = dR["N2m"]
    dvdw["LYS"]["CA"] = dR["CT"]
    dvdw["LYS"]["C"] = dR["C"]
    dvdw["LYS"]["O"] = dR["O"]
    dvdw["LYS"]["CB"] = dR["CT"]
    dvdw["LYS"]["CG"] = dR["CT"]
    dvdw["LYS"]["CD"] = dR["CT"] 
    dvdw["LYS"]["CE"] = dR["CT"]
    dvdw["LYS"]["NZ"] = dR["N3n"] 
    dvdw["LYS"]["H"] = dR["H"]
    dvdw["LYS"]["HA"] = dR["H1"]
    dvdw["LYS"]["HB2"] = dR["HC"]
    dvdw["LYS"]["HB3"] = dR["HC"]
    dvdw["LYS"]["HG2"] = dR["HC"]
    dvdw["LYS"]["HG3"] = dR["HC"]
    dvdw["LYS"]["HD2"] = dR["HC"]
    dvdw["LYS"]["HD3"] = dR["HC"]
    dvdw["LYS"]["HE2"] = dR["HP"]
    dvdw["LYS"]["HE3"] = dR["HP"]
    dvdw["LYS"]["HZ2"] = dR["H"]
    dvdw["LYS"]["HZ1"] = dR["H"]   
    dvdw["LYS"]["HZ3"] = dR["H"]

    dvdw["HIS"] = {}
    dvdw["HIS"]["N"] = dR["N2m"]
    dvdw["HIS"]["CA"] = dR["CT"]
    dvdw["HIS"]["C"] = dR["C"]
    dvdw["HIS"]["O"] = dR["O"]
    dvdw["HIS"]["CB"] = dR["CT"]
    dvdw["HIS"]["CG"] = dR["CA"]
    dvdw["HIS"]["ND1"] = dR["N2m"]
    dvdw["HIS"]["CE1"] = dR["CA"]
    dvdw["HIS"]["NE2"] = dR["N2m"]
    dvdw["HIS"]["CD2"] = dR["CA"]
    dvdw["HIS"]["H"] = dR["H"]
    dvdw["HIS"]["HA"] = dR["H1"]
    dvdw["HIS"]["HB2"] = dR["HC"]
    dvdw["HIS"]["HB3"] = dR["HC"]
    dvdw["HIS"]["HE1"] = dR["H5"]
    dvdw["HIS"]["HE2"] = dR["H"]
    dvdw["HIS"]["HD2"] = dR["H4"]

    dvdw["HID"] = {}
    dvdw["HID"]["N"] = dR["N2m"]
    dvdw["HID"]["CA"] = dR["CT"]
    dvdw["HID"]["C"] = dR["C"]
    dvdw["HID"]["O"] = dR["O"]
    dvdw["HID"]["CB"] = dR["CT"]
    dvdw["HID"]["CG"] = dR["CA"]
    dvdw["HID"]["ND1"] = dR["N2m"]
    dvdw["HID"]["CE1"] = dR["CA"]
    dvdw["HID"]["NE2"] = dR["N2m"]
    dvdw["HID"]["CD2"] = dR["CA"]
    dvdw["HID"]["H"] = dR["H"]
    dvdw["HID"]["HA"] = dR["H1"]
    dvdw["HID"]["HB2"] = dR["HC"]
    dvdw["HID"]["HB3"] = dR["HC"]
    dvdw["HID"]["HE1"] = dR["H5"]
    dvdw["HID"]["HD1"] = dR["H"]
    dvdw["HID"]["HD2"] = dR["H4"]

    dvdw["MET"] = {}
    dvdw["MET"]["N"] = dR["N2m"]
    dvdw["MET"]["CA"] = dR["CT"]
    dvdw["MET"]["C"] = dR["C"]
    dvdw["MET"]["O"] = dR["O"]
    dvdw["MET"]["CB"] = dR["CT"]
    dvdw["MET"]["CG"] = dR["CT"]
    dvdw["MET"]["SD"] = dR["S"]
    dvdw["MET"]["CE"] = dR["CT"]
    dvdw["MET"]["H"] = dR["H"]
    dvdw["MET"]["HA"] = dR["H1"]
    dvdw["MET"]["HB2"] = dR["HC"]
    dvdw["MET"]["HB3"] = dR["HC"]
    dvdw["MET"]["HG2"] = dR["H1"]
    dvdw["MET"]["HG3"] = dR["H1"]
    dvdw["MET"]["HE2"] = dR["H1"]
    dvdw["MET"]["HE3"] = dR["H1"]
    dvdw["MET"]["HE1"] = dR["H1"]

    dvdw["PHE"] = {}
    dvdw["PHE"]["N"] = dR["N2m"]
    dvdw["PHE"]["CA"] = dR["CT"]
    dvdw["PHE"]["C"] = dR["C"]
    dvdw["PHE"]["O"] = dR["O"]
    dvdw["PHE"]["CB"] = dR["CT"]
    dvdw["PHE"]["CG"] = dR["CA"]
    dvdw["PHE"]["CD1"] = dR["CA"]
    dvdw["PHE"]["CD2"] = dR["CA"]
    dvdw["PHE"]["CE1"] = dR["CA"]
    dvdw["PHE"]["CE2"] = dR["CA"]
    dvdw["PHE"]["CZ"] = dR["CA"]
    dvdw["PHE"]["H"] = dR["H"]
    dvdw["PHE"]["HA"] = dR["H1"]
    dvdw["PHE"]["HB2"] = dR["HC"]
    dvdw["PHE"]["HB3"] = dR["HC"]
    dvdw["PHE"]["HD1"] = dR["HA"]
    dvdw["PHE"]["HD2"] = dR["HA"]
    dvdw["PHE"]["HE1"] = dR["HA"]
    dvdw["PHE"]["HE2"] = dR["HA"]
    dvdw["PHE"]["HZ"] = dR["HA"]
   
    dvdw["TYR"] = {}
    dvdw["TYR"]["N"] = dR["N2m"]
    dvdw["TYR"]["CA"] = dR["CT"]
    dvdw["TYR"]["C"] = dR["C"]
    dvdw["TYR"]["O"] = dR["O"]
    dvdw["TYR"]["CB"] = dR["CT"]
    dvdw["TYR"]["CG"] = dR["CA"]
    dvdw["TYR"]["CD1"] =  dR["CA"]
    dvdw["TYR"]["CD2"] =  dR["CA"]
    dvdw["TYR"]["CE1"] = dR["CA"]
    dvdw["TYR"]["CE2"] = dR["CA"]
    dvdw["TYR"]["CZ"] = dR["CA"]
    dvdw["TYR"]["OH"] = dR["OH"]
    dvdw["TYR"]["H"] = dR["H"]
    dvdw["TYR"]["HA"] = dR["H1"]
    dvdw["TYR"]["HB2"] = dR["HC"]
    dvdw["TYR"]["HB3"] = dR["HC"]
    dvdw["TYR"]["HD1"] = dR["HA"]
    dvdw["TYR"]["HD2"] = dR["HA"]
    dvdw["TYR"]["HE1"] = dR["HA"]
    dvdw["TYR"]["HE2"] = dR["HA"]
    dvdw["TYR"]["HH"] = dR["HO"]

    dvdw["TRP"] = {}
    dvdw["TRP"]["N"] = dR["N2m"]
    dvdw["TRP"]["CA"] = dR["CT"]
    dvdw["TRP"]["C"] = dR["C"]
    dvdw["TRP"]["O"] = dR["O"]
    dvdw["TRP"]["CB"] = dR["CT"]
    dvdw["TRP"]["CG"] = dR["CA"]
    dvdw["TRP"]["CD1"] = dR["CA"]
    dvdw["TRP"]["CD2"] = dR["CA"]
    dvdw["TRP"]["NE1"] = dR["N2m"]
    dvdw["TRP"]["CE2"] = dR["CA"]
    dvdw["TRP"]["CE3"] = dR["CA"]
    dvdw["TRP"]["CZ2"] = dR["CA"]
    dvdw["TRP"]["CZ3"] = dR["CA"]
    dvdw["TRP"]["CH2"] = dR["CA"]
    dvdw["TRP"]["H"] = dR["H"]
    dvdw["TRP"]["HA"] = dR["H1"]
    dvdw["TRP"]["HB2"] = dR["HC"]
    dvdw["TRP"]["HB3"] = dR["HC"]
    dvdw["TRP"]["HD1"] = dR["H4"]
    dvdw["TRP"]["HE1"] = dR["H"]
    dvdw["TRP"]["HE3"] = dR["HA"]
    dvdw["TRP"]["HZ2"] = dR["HA"]
    dvdw["TRP"]["HZ3"] = dR["HA"]
    dvdw["TRP"]["HH2"] = dR["HA"]

    if residue_name in dvdw and atom_name in dvdw[residue_name]:
        return dvdw[residue_name][atom_name]
    else:
        # Domylnie zwróæ 0, jeli nie ma informacji o promieniu
        return 0.0

def calculate_clash_score(structure, distance_threshold=0.4):

    n_atoms = 0
    clash_count = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    n_atoms += 1

                    for other_model in structure:
                        for other_chain in other_model:
                            for other_residue in other_chain:
                                for other_atom in other_residue:
                                    if atom != other_atom:

                                        vdw_atom = epsilon_vdw_PDB(atom.get_name(), residue.get_resname())
                                        vdw_other_atom = epsilon_vdw_PDB(other_atom.get_name(), other_residue.get_resname())

                                        if vdw_atom is not None and vdw_other_atom is not None:    

                                            distance = atom - other_atom
                                            vdw_sum = vdw_atom + vdw_other_atom + distance_threshold  

                                            if atom.parent.id == other_atom.parent.id:
                                                continue
                                            # Wykluczamy clashes z wi¹zaniami peptydowymi
                                            elif (other_atom.name == "C" and atom.name == "N") or (other_atom.name == "N" and atom.name == "C"):
                                                continue
                                            # Wykluczamy clashes z mostkami dwusiarczkowymi
                                            elif (other_atom.name == "SG" and atom.name == "SG") and distance > 1.88:
                                                continue

                                            if distance < vdw_sum:
                                                print(atom.get_name(), other_atom.get_name())
                                                print(distance, vdw_sum)
                                                clash_count += 1
    print(clash_count)    
    print(n_atoms)

    # Obliczamy clash score
    clash_score = (clash_count / n_atoms) * 1000
    return clash_score, n_atoms


if __name__ == "__main__":
    pdb_file_path = "C:\\Users\\anna\\Downloads\\101m.pdb"

    structure = PDB.PDBParser(QUIET=True).get_structure("rna_structure", pdb_file_path)

    clash_score, n_atoms = calculate_clash_score(structure, distance_threshold=0.4)

    print(f"Clash Score: {clash_score:.2f}")

