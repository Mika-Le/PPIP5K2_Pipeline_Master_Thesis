from openeye.oechem import *
from kinoml.modeling.OEModeling import read_molecules, read_smiles, prepare_complex, update_residue_identifiers, write_molecules
from kinoml.docking.OEDocking import create_hybrid_receptor, hybrid_docking

import multiprocessing
import sys

#read hits
hit_file = str(sys.argv[1])
ifs = oemolistream("/scratch/mikale/docking/" + hit_file)
# create an empty list
mollist = []

# loop over input
while ifs.IsValid():
    mol = OEGraphMol()
    OEReadMolecule(ifs, mol)
    mollist.append(mol)

print("Number of hits:", len(mollist))

#read protein

design_unit = OEDesignUnit()
OEReadDesignUnit("/scratch/mikale/4nzn_design_unit.oedu", design_unit)

protein, ligand = OEGraphMol(), OEGraphMol()
design_unit.GetProtein(protein)
design_unit.GetLigand(ligand)
hybrid_receptor = create_hybrid_receptor(protein, ligand)

def dock_molecule(molec):
    print("Docking next molecule")
    try:
        docking_poses = hybrid_docking(hybrid_receptor, [molec], num_poses=3)
        return docking_poses
    except:
        print("Docking failed")

#n_cpus = multiprocessing.cpu_count()
#dock molecules

ofs = oemolostream("/scratch/mikale/docking/docking_" + hit_file)
ofs.SetFormat(OEFormat_SDF)
with multiprocessing.Pool(16) as pool:
    dockings = pool.imap(dock_molecule, mollist, chunksize=100)
    if dockings is not None:
        for poses in dockings:
            if poses is not None:
                for pose in poses:
                    OEWriteMolecule(ofs, pose)