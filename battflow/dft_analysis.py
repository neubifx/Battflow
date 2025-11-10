from pathlib import Path
import os
import shutil

from ase.io import read, write

from MDAnalysis.coordinates.XYZ import XYZWriter
from MDAnalysis.topology.guessers import guess_types

from cclib.parser import ccopen

from ase.calculators.orca import ORCA, OrcaProfile

from bson import ObjectId



def create_dft_folders_and_coordinates(m_smiles, mols, a_smiles, ans, c_smiles, cats, ions_pdb, dft_path, pdb_files):

    #create binding energies folder

    binding_energies_folder = Path(dft_path) / "00_binding_energies"
    binding_energies_folder.mkdir(parents=True, exist_ok=True)
    
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

    return  binding_energies_folder, dft_component_folders, packed_molecules, packed_anions, packed_cations
    
    
def extract_solvation_structure_from_md(u, solvation_shell, solute, binding_energies_folder):
    binding_energies_folder = Path(binding_energies_folder)
    binding_energies_folder.mkdir(parents=True, exist_ok=True)

    #  Remove existing subfolders
    for sub in binding_energies_folder.iterdir():
        if sub.is_dir():
            shutil.rmtree(sub)

    n_frames = len(u.trajectory)
    frame_to_use = int(n_frames - 1)

    # Ensure the universe has element symbols for XYZ output
    elem_all = guess_types(u.atoms.names)
    u.add_TopologyAttr("elements", elem_all)

    for i, shell in enumerate(solvation_shell):
        #  Skip keys with 0 count or "fraction"
        shell_dict = {k: v for k, v in shell.items() if k != "fraction" and v != 0}
        print(f"\nGetting shells for: {shell_dict}")

        #  Create output folder
        folder_name = f"{i+1:02d}." + "_".join(f"{k}_{v}" for k, v in shell_dict.items())
        shell_folder = binding_energies_folder / folder_name
        shell_folder.mkdir(parents=True, exist_ok=True)

        #  Get solvation shell indices
        shell_frames = solute.speciation.get_shells(shell_dict)
        print(shell_frames.head())

        first_solute_ix = shell_frames.index.get_level_values("solute_ix")[0]
        print(f"First solute_ix: {first_solute_ix}")

        #try to avoid catching additional residues
        n_closest = sum(shell_dict.values())

        shell_coord = solute.get_shell(solute_index=first_solute_ix, frame=frame_to_use, closest_n_only=n_closest)
        shell_coord.unwrap(compound="residues")

        xyz_path = shell_folder / "shell_coord.xyz"

        #  Write XYZ file
        with XYZWriter(xyz_path, n_atoms=shell_coord.n_atoms) as w:
            w.write(shell_coord)

        print(f" Saved {xyz_path.name} in {shell_folder.name}")

        #  Read content into the shell entry
        with open(xyz_path, "r") as f:
            shell["xyz"] = f.read()

    return solvation_shell
    
    
def extract_orbital_energies(output_file):
    """
    Extract HOMO, LUMO, and HOMO–LUMO gap using cclib from ORCA output.

    Returns a list of dicts: one per spin ('alpha', 'beta').
    """
    data = ccopen(str(output_file)).parse()
    mo_es = data.moenergies
    homos = data.homos
    results = []

    for spin, levels in enumerate(mo_es):
        homo_i = homos[spin]
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
    
    
def calculate_component_dft_energies(dft_component_folders, mols, ans, cats, config, packed_molecules, packed_cations, packed_anions):
    """
    Calculate DFT energies, charges, and orbital energies for each component.

    Returns:
        dict: {label: {'energy': ..., 'charge': ..., 'HOMO': ..., ...}}
    """
    energies_dict = {}

    packed_all = packed_molecules + packed_anions + packed_cations
    labels = mols + ans + cats

    for folder, label, (smiles, _) in zip(dft_component_folders, labels, packed_all):
        folder = Path(folder)
        xyz_path = folder / "residue_coord.xyz"

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
            num_plus = smiles.count('+')
            num_minus = smiles.count('-')
            charge = num_plus - num_minus

        print(f" → Using charge {charge}")

        # Setup ORCA calculator
        calc = ORCA(
            profile=OrcaProfile(command=config["dft_simulations"]["orca_profile"]),
            charge=charge,
            mult=1,
            orcasimpleinput=config["dft_simulations"]["orca_input_block"],
            orcablocks=(
                f'%pal nprocs {config["dft_simulations"]["ncores"]} end '
                f'%geom maxiter 1000 end %scf maxiter 1000 end'
            ),
            directory=str(folder)
        )

        atoms.calc = calc

        try:
            energy = atoms.get_potential_energy()
            print(f" Energy for {label}: {energy:.6f} eV")

            result = {
                "energy": energy,
                "charge": charge,
            }

            orca_out = folder / "orca.out"
            if orca_out.exists():
                orbital_info = extract_orbital_energies(orca_out)
                for spin_result in orbital_info:
                    spin = spin_result["spin"]
                    result["HOMO"] = spin_result["HOMO (eV)"]
                    result["LUMO"] = spin_result["LUMO (eV)"]
                    result["HOMO-LUMO"] = spin_result["Gap (eV)"]
            else:
                print(f" {orca_out.name} not found. HOMO/LUMO not extracted.")

            energies_dict[label] = result

        except Exception as e:
            print(f" Failed for {label}: {e}")
            energies_dict[label] = {"energy": None, "charge": charge}

    return energies_dict
    
    
def calculate_dft_binding_energies(binding_energies_folder, solvation_shell, config, energies_dict):
    """
    Calculates DFT energy and orbital info for each solvation shell.

    Args:
        binding_energies_folder (str or Path): Path to folders with shell_coord.xyz.
        solvation_shell (list of dict): Component counts per shell.
        config (dict): Parsed YAML config.
        energies_dict (dict): Per-species energies and charges.

    Returns:
        list of dicts: Updated solvation_shell entries with 'energy', 'homo', etc.
    """
    binding_energies_folder = Path(binding_energies_folder)
    subfolders = sorted(binding_energies_folder.iterdir())
    updated_shells = []

    for i, (folder, shell_entry) in enumerate(zip(subfolders, solvation_shell)):
        xyz_path = folder / "shell_coord.xyz"

        if not xyz_path.exists():
            print(f" Skipping {folder.name}: shell_coord.xyz not found.")
            shell_entry["energy"] = None
            updated_shells.append(shell_entry)
            continue

        print(f" Running ORCA for shell {i+1} in folder {folder.name}")

        atoms = read(xyz_path)

        #  Total charge starts with +1 from Li+
        total_charge = 1
        for component, count in shell_entry.items():
            if component in ["fraction", "energy", "binding_energy", "xyz"]:
                continue
            if component in energies_dict and energies_dict[component]["charge"] is not None:
                total_charge += count * energies_dict[component]["charge"]
            else:
                print(f" Missing charge info for '{component}' — assuming charge 0")

        print(f" → Total shell charge: {total_charge}")

        #  ORCA calculator setup
        calc = ORCA(
            profile=OrcaProfile(command=config["dft_simulations"]["orca_profile"]),
            charge=total_charge,
            mult=1,
            orcasimpleinput=config["dft_simulations"]["orca_input_block"],
            orcablocks=(
                f'%pal nprocs {config["dft_simulations"]["ncores"]} end '
                f'%geom maxiter 1000 end %scf maxiter 1000 end'
            ),
            directory=str(folder)
        )

        atoms.calc = calc

        try:
            energy = atoms.get_potential_energy()
            print(f" Shell energy: {energy:.6f} eV")
            shell_entry["energy"] = energy

            # Extract HOMO/LUMO/GAP from ORCA output
            orca_out = folder / "orca.out"
            if orca_out.exists():
                orbital_info = extract_orbital_energies(orca_out)
                for spin_result in orbital_info:
                    spin = spin_result["spin"]
                    shell_entry["HOMO"] = spin_result["HOMO (eV)"]
                    shell_entry["LUMO"] = spin_result["LUMO (eV)"]
                    shell_entry["HOMO-LUMO"] = spin_result["Gap (eV)"]
            else:
                print(f" {orca_out.name} not found. HOMO/LUMO not extracted.")

        except Exception as e:
            print(f" ORCA failed in {folder.name}: {e}")
            shell_entry["energy"] = None

        updated_shells.append(shell_entry)

    return updated_shells
    
def compute_binding_energies(updated_shells, energies_dict, config):
    """
    Compute binding energy for each solvation shell and update each entry in-place.
    
    Args:
        updated_shells (list of dict): Each dict contains component counts, shell DFT energy, etc.
        energies_dict (dict): Per-component info, e.g. {'mp': {'energy': -7300, 'charge': 0, ...}}
        config (dict): DFT config dict with Li energy under ["dft_simulations"]["li_energy"]

    Returns:
        list of dict: updated_shells with new key 'binding_energy' (in eV)
    """
    li_energy = config["dft_simulations"]["li_energy"]

    # Keys that are not molecular components
    exclude_keys = {"energy", "fraction", "binding_energy", "HOMO", "LUMO", "HOMO-LUMO", "xyz"}

    for entry in updated_shells:
        shell_energy = entry.get("energy")
        if shell_energy is None:
            entry["binding_energy"] = None
            continue

        # Start with E(shell) - E(Li)
        binding_energy = shell_energy - li_energy

        for component, count in entry.items():
            if component in exclude_keys:
                continue

            if component not in energies_dict:
                print(f" Warning: {component} not found in energies_dict, skipping.")
                continue

            component_energy = energies_dict[component].get("energy")
            if component_energy is None:
                print(f" Missing energy for {component}, skipping.")
                continue

            binding_energy -= count * component_energy

        entry["binding_energy"] = binding_energy

    return updated_shells
    
def upload_dft_calculated_data(doc_id, updated_shells, energies_dict, collection):
    collection.update_one(
        {"_id": ObjectId(doc_id)},
        {"$set": {
            "simulation_data.solvation_statistics": updated_shells,
            "simulation_data.dft_energies": energies_dict
        }}
    )

    
