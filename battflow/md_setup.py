from pathlib import Path
import os
import shutil

import acpype
from acpype.topol import ACTopol, MolTopol

import MDAnalysis as mda
import mdapackmol


def prepare_simulation_paths(doc_id):
    """
    Prepare working and setup paths for a simulation based on a document ID.

    Args:
        doc_id (ObjectId or str): Unique identifier of the document being processed.

    Returns:
        tuple: (work_path, setup_path, ff_path, pack_path, md_path, dft_path)
            - work_path (Path): Current working directory.
            - setup_path (Path): Directory path for the current simulation setup.
            - ff_path (Path): Force fields path.
            - pack_path (Path): Electrolyte structure path.
            - md_path (Path): Molecular dynamics run path.
            - dft_path (Path): DFT simulations path.
    """
    work_path = Path.cwd()
    setup_path = work_path / "working_md_dir" / str(doc_id)
    ff_path = setup_path / "00.force_fields"
    pack_path = setup_path / "01.electrolyte_structure"
    md_path = setup_path / "02.md_run"
    dft_path = setup_path / "03.dft_components"

    # Create main setup directory
    setup_path.mkdir(parents=True, exist_ok=True)

    # Create subfolders
    folders = [ff_path, pack_path, md_path, dft_path]
    for folder in folders:
        if folder.exists():
            print(f"The folder {folder.name} already exists!")
        else:
            print(f"Creating {folder.name} ...")
            folder.mkdir()

    return work_path, setup_path, ff_path, pack_path, md_path, dft_path
    

def names_smiles_molarity_setup(doc):
    """
    Returns the names and SMILES strings of molecules and anions 
    from a given electrolyte document.

    Args:
        doc (dict): A single document from the electrolyte database, containing
            molecule and anion names and their corresponding SMILES strings.

    Returns:
         mols (list): Names of molecules
         ans (list): Names of anions
         m_smiles (list): SMILES strings of molecules
         a_smiles (list): SMILES strings of anions
    """
            
    #electrolyte composition
    mols = doc["components"]["molecules"]
    ans = doc["components"]["anions"]
    cats = doc["components"]["cations"]
    ions = doc["components"]["ions"]
    
    #smiles
    m_smiles = doc["smiles"]["molecules"]
    a_smiles = doc["smiles"]["anions"]
    c_smiles = doc["smiles"]["cations"]

    #concentrations
    m_conc = doc["concentrations"]["molecules"]
    a_conc = doc["concentrations"]["anions"]
    c_conc = doc["concentrations"]["cations"]
    i_conc = doc["concentrations"]["ions"]
    
    return mols, ans, cats, ions, m_smiles, a_smiles, c_smiles, m_conc, a_conc, c_conc, i_conc
    
    
def prepare_molecule_topologies(work_path, ff_path, mols, m_smiles):
    """
    Generate topology files for each molecule using ACPYPE.

    Args:
        work_path (str): Path to return to after setup.
        setup_path (str): Path where component folders will be created.
        mols (list): List of molecule names.
        m_smiles (list): List of SMILES strings corresponding to mols.
    """
    if not mols:
        print("No molecules found in the electrolyte. Skipping topology generation") #Only if mols is empty
        return
    
    for i, m in enumerate(m_smiles):
        folder_name = f"{i+1:02d}.{mols[i]}_{m}"
        component_folder = Path(ff_path) / folder_name

        if component_folder.exists():
            print(f"{folder_name} folder already exists")
        else:
            component_folder.mkdir(parents = True)

        # Enter in the folder respective to the component
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

        # Edit resname in .mol2 file in order starting from 001
        mol2_file = Path(component_folder) / f"{mols[i]}.mol2"
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
        

def prepare_anion_topologies(work_path, ff_path, ans, a_smiles, m_smiles):
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

    if not ans:
        print("No anions found in the electrolyte. Skipping topology generation") #Only if ans is empty
        return
    
    offset = len(m_smiles) # enumerated after molecules topologies

    for i, a in enumerate(a_smiles):
        index = i + 1 + offset
        folder_name = f"{index:02d}.{ans[i]}_{a}"
        component_folder = Path(ff_path) / folder_name

        if component_folder.exists():
            print(f"{ans[i]}_{a} folder already exists")
        else:
            component_folder.mkdir(parents = True)

        # Enter in the folder respective to the component
        os.chdir(component_folder)

        # Count charges from SMILES
        charge = a.count("-")
        multiplicity = 1 if charge == 2 else 2

        # Generate anion topology
        molecule = ACTopol(
            a,
            chargeType="bcc",
            basename=ans[i],
            chargeVal=-charge,
            multiplicity=multiplicity,
            verbose=True
        )

        # Modify .mol2 file
        mol2_file = Path(component_folder) / f"{ans[i]}.mol2"
        with open(mol2_file, "r") as f:
            lines = f.readlines()

        # Rename file to follow the resname order, considering the offset
        lines = [line.replace("UNL1", f"{index:03d}") for line in lines] 

        with open(mol2_file, "w") as f:
            f.writelines(lines)

        # Generate topology files
        molecule.createACTopol()
        molecule.createMolTopol()

        # Return to the working directory
        os.chdir(work_path)
        

def prepare_cation_topologies(work_path, ff_path, cats, c_smiles, a_smiles, m_smiles):
    """
    Generate topology files for each cation using ACPYPE.

    Args:
        work_path (str): The base working directory to return to after processing.
        setup_path (str): Path where anion folders will be created.
        cats (list): List of cations names.
        c_smiles (list): Corresponding list of SMILES strings for cations.
        a_smiles (list): Corresponding list of SMILES strings for anions.
        m_smiles (list): List of molecule SMILES used to offset index count.

    Returns:
        None
    """

    if not cats:
        print("No cations found in the electrolyte. Skipping topology generation") #Only if ans is empty
        return
    
    offset = len(m_smiles + a_smiles) # enumerated after molecules topologies

    for i, c in enumerate(c_smiles):
        index = i + 1 + offset
        folder_name = f"{index:02d}.{cats[i]}_{c}"
        component_folder = Path(ff_path) / folder_name

        if component_folder.exists():
            print(f"{cats[i]}_{c} folder already exists")
        else:
            component_folder.mkdir(parents = True)

        # Enter in the folder respective to the component
        os.chdir(component_folder)

        # Count charges from SMILES
        charge = c.count("+")
        multiplicity = 1 if charge ==2 else 2

        # Generate anion topology
        molecule = ACTopol(
            c,
            chargeType="bcc",
            basename=cats[i],
            chargeVal=charge,
            multiplicity=multiplicity,
            verbose=True
        )

        # Modify .mol2 file
        mol2_file = Path(component_folder) / f"{cats[i]}.mol2"
        with open(mol2_file, "r") as f:
            lines = f.readlines()

        # Rename file to follow the resname order, considering the offset
        lines = [line.replace("UNL1", f"{index:03d}") for line in lines] 

        with open(mol2_file, "w") as f:
            f.writelines(lines)

        # Generate topology files
        molecule.createACTopol()
        molecule.createMolTopol()

        # Return to the working directory
        os.chdir(work_path)


def process_ion_topologies(BASE_DIR, config, ions, pack_path, md_path):
    """
    Identify and copy ion topologies files to their respective folders

    Args:
        config (dict): Dictionary containing the load .yaml file
        ions (list): Names of ions
        pack_path (Path): Electrolyte structure path.

    Returns:
        ions_itp_file (pathlib.Path): Topology file for ions
        topol_main_file (pathlib.Path): Main topol file containing ions topologies
        ions_pdb (list): List containing the pathlib.Path for the pdb files for ions 
        
    """
    
    
    topol_main_file = BASE_DIR / config["md_simulations"]["topol_main"]
    li_pdb = BASE_DIR / config["md_simulations"]["li_pdb"]
    
    #define dict to search for ion pdb files. Add other ions later.    
    ion_files = {
        "li" : li_pdb
    }

    ions_pdb = []
    for ion, pdb_file in ion_files.items():
        if ion in ions:
            shutil.copy(pdb_file, pack_path) 
            ions_pdb.append(pdb_file)
            
            #Separating the ions itp file and processing in this loop
            ions_itp_file = BASE_DIR / config["md_simulations"][f"{ion}_itp"]
            shutil.copy(ions_itp_file, md_path) 
            
    shutil.copy(topol_main_file, md_path) 
    
    return ions_itp_file, topol_main_file, ions_pdb
    

def process_all_topologies(m_smiles, mols, a_smiles, ans, c_smiles, cats, ions_pdb, ff_path, pack_path, md_path):
    """
    Function to move itp, top and pdb files to their respective folders for simulation setup

    Args:
        mols (list): Names of molecules
        ans (list): Names of anions
        m_smiles (list): SMILES strings of molecules
        a_smiles (list): SMILES strings of anions
    """
    
    packed_molecules = [(m_smiles[i], mols[i]) for i in range(len(m_smiles))] 
    packed_anions = [(a_smiles[i], ans[i]) for i in range(len(a_smiles))] 
    packed_cations = [(c_smiles[i], cats[i]) for i in range(len(c_smiles))]
    
    combined_packed = packed_molecules + packed_anions + packed_cations

    mols_ans_cats_files = []
    itp_files = []
    top_files = []
    for i, (smiles, label) in enumerate(combined_packed):
        
        folder_name = f"{i+1:02d}.{label}_{smiles}"
        acpype_folder = Path(ff_path) / folder_name / f"{label}.acpype"
        print(f"Processing files in {os.path.basename(acpype_folder)}")
        
        #defining relevant files
        pdb_file = acpype_folder / f"{label}_NEW.pdb"
        itp_file = acpype_folder / f"{label}_GMX.itp"
        top_file = acpype_folder / f"{label}_GMX.top"  

        #Append file paths into lists to preserve ordering
        mols_ans_cats_files.append(pdb_file)    
        itp_files.append(itp_file)
        top_files.append(top_file)    
        
        #copying files to desirable folders
        shutil.copy(pdb_file, pack_path)
        shutil.copy(itp_file, md_path)
        #shutil.copy(top_file, md_path) #not needed

    pdb_files = mols_ans_cats_files + ions_pdb

    return pdb_files, itp_files, top_files
    
def number_of_molecules(m_conc, a_conc, c_conc, i_conc):
    """
    Define the number of molecules for each component based on their concentration

    Args:
        m_conc (list) : Concentration of molecules
        a_conc (list) : Concentration of anions
        c_conc (list) : Concentration of cations
        i_conc (list) : Concentration of ions

    Returns:
        n_mol_box (list): List containing the number of molecules for each component
    """

    #box settings 
    a_side = 50 # fixed 50 Angstroms box
    box_vol = (a_side**3) * 1e-27 # box volume in L
    avo = 6.022e23
    
    combined_concs = m_conc + a_conc + c_conc + i_conc
    
    n_mols_box = [] #list containing the number of molecules for each component
    
    for conc in combined_concs:
        conc_value = float(conc.replace("M", "")) #removes concentration from the string and turn it into a value
        n_mol = round(conc_value * box_vol * avo)
        n_mols_box.append(n_mol)  
    
    return a_side, n_mols_box
    
    
def packmol_build(work_path, pack_path, md_path, pdb_files, a_side, n_mols_box, ions):
    """
    Create the electrolyte structure using MDAPackmol

    Args:
        work_path (Pathlib): path for the work dir
        pack_path (Pathlib): path for the electrolyte structure folder
        pdb_files (list): list with paths for all pdb files
        n_mol_box (list): List containing the number of molecules for each component
        
    Returns:
        system (MDAUniverse): MDA Universe object containing electrolyte structure
        packmol_file (Pathlib) : Path for the .gro file created
    
    """

    temp_folder = ("/home/neubijr/test_packmol") #used to avoid problems during script generation  
    os.chdir(temp_folder) #only in script testing due to permission problems
    #os.chdir(pack_path) #use in final release

    packed_concentrations = [(pdb_files[i], n_mols_box[i]) for i in range(len(pdb_files))]
    
    packmol_system = []
    #create electrolyte box based on MDAnalysis Universe
    for (path, number) in packed_concentrations:
        comp = mdapackmol.PackmolStructure(mda.Universe(path), number, 
                                    instructions=[f"inside box 0. 0. 0. {a_side}. {a_side}. {a_side}."])
        packmol_system.append(comp)

    #build
    system = mdapackmol.packmol(packmol_system)
    system.dimensions = [a_side, a_side, a_side, 90, 90, 90]

    #new loop to define resnames that get lost after packing

    #Extremely rigid as it is mandatory to keep name of ions as li.pdb, na.pdb etc. Needs to be changed later.
    
    res_index = 0
    for i, (pdb_file, amount) in enumerate(packed_concentrations): 
        
        filename = Path(pdb_file).stem #Use .stem from Path to get the file name without extension
        for j in range(amount):   
            
            if str(filename) in ions:
                system.residues[res_index].resname = str(filename)
            else:
                system.residues[res_index].resname = f"00{i+1}" 

            res_index += 1
        

    #save electrolyte gro file
    packmol_file = Path(pack_path) / "electrolyte.gro"
    system.atoms.write(packmol_file)

    #copy to md_run folder
    shutil.copy(packmol_file, md_path) 
    
    os.chdir(work_path)

    return system, packmol_file

