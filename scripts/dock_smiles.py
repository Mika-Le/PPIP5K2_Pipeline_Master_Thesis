from openeye import oechem

from kinoml.modeling.OEModeling import read_molecules, read_smiles, prepare_complex, update_residue_identifiers, write_molecules
from kinoml.docking.OEDocking import create_hybrid_receptor, hybrid_docking

#read protein and ligand
design_unit = oechem.OEDesignUnit()
oechem.OEReadDesignUnit("../data/4nzn_design_unit.oedu", design_unit)

protein, ligand = oechem.OEGraphMol(), oechem.OEGraphMol()
design_unit.GetProtein(protein)
design_unit.GetLigand(ligand)

#SMILES and names of ligands for docking
print("Preparing small molecules for docking ...")
molecules = [
    {"name": "105229",
     "smiles": "O=C(NC=N1)C2=C1C=C(C(NC3=NC=CS3)=O)C=C2"},
]

oemols = []
for molecule in molecules:
    oemol = read_smiles(molecule["smiles"])
    oemol.SetTitle(molecule["name"])
    oemols.append(oemol)

#docking

print("Docking small molecules ...")
hybrid_receptor = create_hybrid_receptor(protein, ligand)
docking_poses = hybrid_docking(hybrid_receptor, oemols, num_poses=20)
write_molecules(docking_poses, "../data/docking_poses.sdf")