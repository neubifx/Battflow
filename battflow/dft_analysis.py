from pathlib import Path
import os
import shutil

from ase.io import read, write

from MDAnalysis.coordinates.XYZ import XYZWriter
from MDAnalysis.topology.guessers import guess_types

from cclib.parser import ccopen

from ase.calculators.orca import ORCA, OrcaProfile

from bson import ObjectId

from battflow.utils import normalize_id

from io import StringIO


def create_dft_folders_and_coordinates(m_smiles, mols, a_smiles, ans, c_smiles, cats, ions_pdb, dft_path, pdb_files):
   
    packed_molecules = [(m_smiles[i], mols[i]) for i in range(len(m_smiles))] 
    packed_anions = [(a_smiles[i], ans[i]) for i in range(len(a_smiles))] 
    packed_cations = [(c_smiles[i], cats[i]) for i in range(len(c_smiles))]

    combined_packed = packed_molecules + packed_anions + packed_cations

    dft_component_folders = []
    for i, (smiles, label) in enumerate(combined_packed):
        
        folder_name = f"{i+1:02d}.{label}_{smiles}"
        dft_component_folder = Path(dft_path) / folder_name 

        # Create folder
        dft_component_folder.mkdir(parents=True, exist_ok=True)

        print(f"Creating files in {dft_component_folder.name}")

        dft_component_folders.append(dft_component_folder)

        pdb_path = Path(pdb_files[i])
        xyz_path = dft_component_folder / "residue_coord.xyz"

        atoms = read(pdb_path)          # Read PDB
        write(xyz_path, atoms)          # Write XYZ  

        print(f"Saved {xyz_path.name}")     

    return  dft_component_folders, packed_molecules, packed_anions, packed_cations
    
    
def extract_solvation_structure_from_md(u, solvation_shell, solute):

    n_frames = len(u.trajectory)
    frame_to_use = int(n_frames - 1)

    u.trajectory[frame_to_use]

    # Ensure the universe has element symbols for XYZ output
    elem_all = guess_types(u.atoms.names)
    u.add_TopologyAttr("elements", elem_all)

    for i, shell in enumerate(solvation_shell):
        #  Skip keys with 0 count or "fraction" or "xyz"
        shell_dict = {k: v for k, v in shell.items() if k not in {"fraction", "xyz"} and v != 0}
        print(f"\nGetting shells for: {shell_dict}")

        #  Get solvation shell indices
        shell_frames = solute.speciation.get_shells(shell_dict)
        print(shell_frames.head())

        first_solute_ix = shell_frames.index.get_level_values("solute_ix")[0]
        print(f"First solute_ix: {first_solute_ix}")

        #try to avoid catching additional residues
        n_closest = sum(shell_dict.values())

        shell_coord = solute.get_shell(solute_index=first_solute_ix, frame=frame_to_use, closest_n_only=n_closest)
        shell_coord.unwrap(compound="residues")

        # Writing XYZ file in memory
        temp_xyz = StringIO()
        with XYZWriter(temp_xyz, n_atoms=shell_coord.n_atoms) as w:
            w.write(shell_coord)
            xyz_str = temp_xyz.getvalue() 
            
        #save extracted xyz into solvation_shell for future dft calculations
        solvation_shell[i]["xyz"] = xyz_str

    return solvation_shell
    
    
def extract_orbital_energies(output_file):
    data = ccopen(str(output_file)).parse()

    mo_es = data.moenergies
    homos = data.homos
    results = []

    for spin, levels in enumerate(mo_es):

        # cclib may give only one HOMO index
        if spin < len(homos):
            homo_i = homos[spin]
        else:
            homo_i = homos[0]   # fallback to alpha HOMO

        if homo_i is None or homo_i >= len(levels):
            homo_e = None
            lumo_e = None
            gap = None
        else:
            homo_e = levels[homo_i]
            lumo_e = levels[homo_i + 1] if homo_i + 1 < len(levels) else None
            gap = (lumo_e - homo_e) if lumo_e is not None else None

        results.append({
            "spin": "alpha" if spin == 0 else "beta",
            "HOMO (eV)": homo_e,
            "LUMO (eV)": lumo_e,
            "Gap (eV)": gap
        })

    return results

    
    
def calculate_component_dft_energies(
    dft_component_folders,
    mols,
    ans,
    cats,
    config,
    packed_molecules,
    packed_cations,
    packed_anions,
):
    """
    Calculate DFT energies, charges, and (optionally) orbital energies
    for each component.

    Returns:
        dict: {label: {'energy': ..., 'charge': ..., 'HOMO': ..., 'LUMO': ..., ...}}
    """
    energies_dict = {}

    packed_all = packed_molecules + packed_anions + packed_cations
    labels = mols + ans + cats

    for folder, label, (smiles, _) in zip(dft_component_folders, labels, packed_all):
        folder = Path(folder)
        xyz_path = folder / "residue_coord.xyz"

        # Sanity check
        if not xyz_path.exists():
            print(f" Skipping {label}: {xyz_path.name} not found.")
            energies_dict[label] = {"energy": None, "charge": None}
            continue

        print(f" Running ORCA for: {label} in folder {folder.name}")

        atoms = read(xyz_path)

        # Determine charge
        if (smiles, label) in packed_molecules:
            charge = 0
        else:
            num_plus = smiles.count("+")
            num_minus = smiles.count("-")
            charge = num_plus - num_minus

        print(f" → Using charge {charge}")

        # Setup ORCA calculator
        calc = ORCA(
            profile=OrcaProfile(command=config["dft_simulations"]["orca_profile"]),
            charge=charge,
            mult=1,
            orcasimpleinput=config["dft_simulations"]["orca_input_block"],
            orcablocks=(
                f"%pal nprocs {config['dft_simulations']['ncores']} end "
                f"%geom maxiter 1000 end "
                f"%scf maxiter 1000 end"
            ),
            directory=str(folder),
        )

        atoms.calc = calc

        # ENERGY CALCULATION. It will kill the job if it fails.
        try:
            energy = atoms.get_potential_energy()
            print(f" Energy for {label}: {energy:.6f} eV")
        except Exception as e:
            print(f" Energy calculation failed for {label}: {e}")
            energies_dict[label] = {"energy": None, "charge": charge}
            continue

        # Store results
        result = {
            "energy": energy,
            "charge": charge,
        }


        # ORBITAL PARSING
        orca_out = folder / "orca.out"
        if orca_out.exists():
            try:
                orbital_info = extract_orbital_energies(orca_out)
                for spin_result in orbital_info:
                    # overwriting
                    result["HOMO"] = spin_result.get("HOMO (eV)")
                    result["LUMO"] = spin_result.get("LUMO (eV)")
                    result["HOMO-LUMO"] = spin_result.get("Gap (eV)")
            except Exception as e:
                print(f" Warning: orbital parsing failed for {label}: {e}")
        else:
            print(f" Warning: {orca_out.name} not found. HOMO/LUMO not extracted.")

        energies_dict[label] = result

    return energies_dict
    
    
def upload_dft_calculated_data(doc_id, solvation_shell, energies_dict, collection):

    safe_id = normalize_id(doc_id)

    collection.update_one(
        {"_id": safe_id},
        {"$set": {
            "simulation_data.solvation_statistics": solvation_shell,
            "simulation_data.dft_energies": energies_dict
        }}
    )

    

