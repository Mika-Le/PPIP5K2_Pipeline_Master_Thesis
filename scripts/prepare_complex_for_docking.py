# KinoML (https://github.com/openkinome/kinoml/commit/d56c3d6bf11614a3ad27d623446ad6586a22341c)
# OpenEye toolkit 2020.2.0
from openeye import oechem

from kinoml.modeling.OEModeling import read_molecules, read_smiles, prepare_complex, update_residue_identifiers, write_molecules
from kinoml.docking.OEDocking import create_hybrid_receptor, hybrid_docking

print("Preparing complex ...")
structure = read_molecules("../data/4nzn.pdb")[0]
for atom in structure.GetAtoms():
    residue_name = oechem.OEAtomGetResidue(atom).GetName().strip()
    if residue_name == "2OU":  # delete second ligand
        structure.DeleteAtom(atom)
    elif residue_name == "ANP":  # only keep adenine core
        atom_name = atom.GetName().strip()
        adenine_atom_names = ["C2", "C4", "C5", "C6", "C8", "N1", "N3", "N6", "N7", "N9"]
        if atom_name not in adenine_atom_names:
            structure.DeleteAtom(atom)
design_unit = prepare_complex(structure, ligand_name='ANP')
oechem.OEWriteDesignUnit("../data/4nzn_design_unit.oedu", design_unit)

print("Saving results ...")
solvated_protein = oechem.OEGraphMol()
design_unit.GetComponents(solvated_protein, oechem.OEDesignUnitComponents_Protein | oechem.OEDesignUnitComponents_Solvent)
solvated_protein = update_residue_identifiers(solvated_protein)
write_molecules([solvated_protein], "../data/4nzn_solvated.pdb")

print("Finished!")


