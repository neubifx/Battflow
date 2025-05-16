import os

import acpype
from acpype.topol import ACTopol, MolTopol

def smiles_setup(doc):
    """
    Returns the names and SMILES strings of molecules and anions 
    from a given electrolyte document.

    Args:
        doc (dict): A single document from the electrolyte database, containing
            molecule and anion names and their corresponding SMILES strings.

    Returns:
        tuple: A tuple of four lists:
            - mols (list): Names of molecules
            - ans (list): Names of anions
            - m_smiles (list): SMILES strings of molecules
            - a_smiles (list): SMILES strings of anions
    """
            
    #electrolyte composition
    mols = doc["components"]["molecules"]
    ans = doc["components"]["anions"]
    
    #smiles
    m_smiles = doc["smiles"]["molecules"]
    a_smiles = doc["smiles"]["anions"]
    
    return mols, ans, m_smiles, a_smiles
    

def prepare_simulation_paths(doc_id):
    """
    Prepare working and setup paths for a simulation based on a document ID.

    Args:
        doc_id (ObjectId or str): Unique identifier of the document being processed.

    Returns:
        tuple: (work_path, setup_path)
            - work_path (str): Current working directory.
            - setup_path (str): Directory path for the current simulation setup.
    """
    work_path = os.getcwd()
    os.chdir(work_path)  # optional; might be redundant if you're already in work_path

    setup_path = os.path.join(work_path, "working_md_dir", str(doc_id))

    return work_path, setup_path


def prepare_molecule_topologies(work_path, setup_path, mols, m_smiles):
    """
    Generate topology files for each molecule using ACPYPE.

    Args:
        work_path (str): Path to return to after setup.
        setup_path (str): Path where component folders will be created.
        mols (list): List of molecule names.
        m_smiles (list): List of SMILES strings corresponding to mols.
    """
    for i, m in enumerate(m_smiles):
        folder_name = f"{i+1:02d}.{mols[i]}_{m}"
        component_folder = os.path.join(setup_path, "00.force_fields", folder_name)

        if os.path.exists(component_folder):
            print(f"{folder_name} folder already exists")
        else:
            os.makedirs(component_folder)

        os.chdir(component_folder)

        # Generate topology
        molecule = ACTopol(
            m,
            chargeType="bcc",
            basename=mols[i],
            chargeVal=0,
            multiplicity=1,
            verbose=True
        )

        # Edit resname in .mol2 file
        mol2_file = os.path.join(component_folder, f"{mols[i]}.mol2")
        with open(mol2_file, "r") as f:
            lines = f.readlines()

        lines = [line.replace("UNL1", f"{i+1:03d}") for line in lines]

        with open(mol2_file, "w") as f:
            f.writelines(lines)

        # Generate remaining files
        molecule.createACTopol()
        molecule.createMolTopol()

        # Return to work path
        os.chdir(work_path)
        

def prepare_anion_topologies(work_path, setup_path, ans, a_smiles, m_smiles):
    """
    Generate topology files for each anion using ACPYPE.

    Args:
        work_path (str): The base working directory to return to after processing.
        setup_path (str): Path where anion folders will be created.
        ans (list): List of anion names.
        a_smiles (list): Corresponding list of SMILES strings.
        m_smiles (list): List of molecule SMILES used to offset index count.

    Returns:
        None
    """
    offset = len(m_smiles)

    for i, a in enumerate(a_smiles):
        index = i + 1 + offset
        folder_name = f"{index:02d}.{ans[i]}_{a}"
        component_folder = os.path.join(setup_path,"00.force_fields", folder_name)

        if os.path.exists(component_folder):
            print(f"{ans[i]}_{a} folder already exists")
        else:
            os.makedirs(component_folder)

        os.chdir(component_folder)

        # Count formal charge from SMILES
        charge = a.count("-")
        multiplicity = 1 if charge in (0, 2) else 2

        # Generate molecule topology
        molecule = ACTopol(
            a,
            chargeType="bcc",
            basename=ans[i],
            chargeVal=charge,
            multiplicity=multiplicity,
            verbose=True
        )

        # Modify .mol2 file
        mol2_file = os.path.join(component_folder, f"{ans[i]}.mol2")
        with open(mol2_file, "r") as f:
            lines = f.readlines()

        lines = [line.replace("UNL1", f"{index:03d}") for line in lines]

        with open(mol2_file, "w") as f:
            f.writelines(lines)

        # Generate topology files
        molecule.createACTopol()
        molecule.createMolTopol()

        # Return to the working directory
        os.chdir(work_path)
