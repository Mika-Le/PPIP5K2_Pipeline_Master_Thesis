"""
Inspried from https://github.com/openkinome/kinoml/commit/550e902037a62c120f71a83dd04f95e2ad408666
"""
import multiprocessing
import pathlib
from typing import  List

from openeye import oechem


def read_smiles(smiles: str) -> oechem.OEGraphMol:
    """
    Read molecule from a smiles string.
    Parameters
    ----------
    smiles: str
        Smiles string.
    Returns
    -------
    molecule: oechem.OEGraphMol
        A molecule as OpenEye molecules.
    """
    ims = oechem.oemolistream()
    ims.SetFormat(oechem.OEFormat_SMI)
    ims.openstring(smiles)

    molecules = []
    for molecule in ims.GetOEMols():
        molecules.append(oechem.OEGraphMol(molecule))

    return molecules[0]


def generate_tautomers(molecule: oechem.OEGraphMol) -> List[oechem.OEGraphMol]:
    """
    Generate reasonable tautomers of a given molecule.
    Parameters
    ----------
    molecule: oechem.OEGraphMol
        An OpenEye molecule.
    Returns
    -------
    tautomers: list of oechem.OEGraphMol
        A list of OpenEye molecules holding the tautomers.
    """
    from openeye import oechem, oequacpac

    tautomer_options = oequacpac.OETautomerOptions()
    tautomer_options.SetMaxTautomersGenerated(4096)
    tautomer_options.SetMaxTautomersToReturn(16)
    tautomer_options.SetCarbonHybridization(True)
    tautomer_options.SetMaxZoneSize(50)
    tautomer_options.SetApplyWarts(True)
    pKa_norm = True
    tautomers = [
        oechem.OEGraphMol(tautomer)
        for tautomer in oequacpac.OEGetReasonableTautomers(
            molecule, tautomer_options, pKa_norm
        )
    ]
    return tautomers


def generate_enantiomers(
        molecule: oechem.OEGraphMol,
        max_centers: int = 12,
        force_flip: bool = False,
        enumerate_nitrogens: bool = False, #No nitrogen enantiomers
) -> List[oechem.OEGraphMol]:
    """
    Generate enantiomers of a given molecule.
    Parameters
    ----------
    molecule: oechem.OEGraphMol
        An OpenEye molecule.
    max_centers: int
        The maximal number of stereo centers to enumerate.
    force_flip: bool
        If specified stereo centers should be ignored.
    enumerate_nitrogens: bool
        If nitrogens with invertible pyramidal geometry should be enumerated.
    Returns
    -------
    enantiomers: list of oechem.OEGraphMol
        A list of OpenEye molecules holding the enantiomers.
    """
    from openeye import oechem, oeomega

    enantiomers = [
        oechem.OEGraphMol(enantiomer)
        for enantiomer in oeomega.OEFlipper(
            molecule, max_centers, force_flip, enumerate_nitrogens, True
        )
    ]
    return enantiomers


def enumerate_smiles(counter, file_path):
    """
    Enumerate undefinded stereo centers and reasonable protonation states, and
    save to file.
    """
    writing_buffer = []
    with open(file_path) as rf:
        #with open(f"/scratch/mikale/screening/enumerated_smiles/enumerated_smiles_{counter}.smi", "w") as wf:
        with open(f"enumerated_smiles_{counter}.smi", "w") as wf:
            for i, line in enumerate(rf.readlines()[1:]):
                if i % 1000 == 0:
                    print(i)
                smiles = line.split("\t")[0]
                name = line.split("\t")[2]
                try:
                    molecule = read_smiles(smiles)
                    oechem.OEDeleteEverythingExceptTheFirstLargestComponent(molecule)
                    molecule.SetTitle(name)
                    tautomers = generate_tautomers(molecule)
                    enantiomers = [generate_enantiomers(tautomer, enumerate_nitrogens=False) for tautomer in tautomers]
                    for tautomers in enantiomers:
                        for tautomer in tautomers:
                            writing_buffer.append(f"{oechem.OEMolToSmiles(tautomer)} {tautomer.GetTitle()}")
                except IndexError:
                    writing_buffer.append(f"{smiles} {name}")

            wf.write("\n".join(writing_buffer))

file_paths = pathlib.Path("/scratch/mikale/screening/molport").glob("*.txt")
n_cpus = multiprocessing.cpu_count()
jobs = [(counter, file_path) for counter, file_path in enumerate(file_paths)]
pool = multiprocessing.Pool(n_cpus)  
pool.starmap(enumerate_smiles, jobs)

#enumerate_smiles(1,pathlib.Path("../testsmiles.txt"))
