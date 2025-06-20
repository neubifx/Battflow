import os

from battflow.database import db_connection, scan_properties_collection, config_path
from battflow.md_setup import (
    # Config and structure setup
    names_smiles_molarity_setup,
    prepare_simulation_paths,

    # Molecule and ion topologies
    prepare_molecule_topologies,
    prepare_anion_topologies,
    prepare_cation_topologies,
    process_ion_topologies,
    process_all_topologies,

    # Packing
    number_of_molecules,
    packmol_build,
)

from battflow.topol_tools import prepare_topol
from battflow.md_run import md_simulation_run


def main():

    parser = argparse.ArgumentParser(description="Run Battflow workflow.")
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to alternative YAML config file (default: config/default.yaml)",
    )
    args = parser.parse_args()
    
    BASE_DIR, config = config_path()
    
    print("Connecting to DB ...")
    db, collection = db_connection(config)
    
    print("Scanning collections for missing properties ...")
    flag, doc_id = scan_properties_collection(collection)
    
    if flag:
        
        print(f"Found missing property in document ID {doc_id}!")
        
        doc = collection.find_one({"_id" : doc_id})
        
        print("Building up electrolyte structure ...")
        
        mols, ans, cats, ions, m_smiles, a_smiles, c_smiles, m_conc, a_conc, c_conc, i_conc = names_smiles_molarity_setup(doc)
        
        print("Molecules:", mols)
        print("Anions:", ans)
        print("Cations:", cats)
        print("Ions:", ions)
        print("Setting up simulation folders...")
        
        work_path, setup_path, ff_path, pack_path, md_path, md_em_path, md_eq_path, md_prod_path, dft_path = prepare_simulation_paths(doc_id)
        
        print("Done!\n")
        print("#################################")
        print("\nPreparing molecules topologies ...")
        print("\n#################################\n")
        
        prepare_molecule_topologies(work_path, ff_path, mols, m_smiles)
        
        print("\nDone!\n")
        print("#################################")
        print("\nPreparing anions topologies ...")
        print("\n#################################\n")
        
        prepare_anion_topologies(work_path, ff_path, ans, a_smiles, m_smiles) 

        print("\nDone!\n")
        print("#################################")
        print("\nPreparing cations topologies ...")
        print("\n#################################\n")

        prepare_cation_topologies(work_path, ff_path, cats, c_smiles, a_smiles, m_smiles)

        print("\nDone!\n")
        
        print("Processing topologies files ...")        
        ions_itp_file, topol_main_file, ions_pdb = process_ion_topologies(BASE_DIR, config, ions, pack_path, md_em_path, md_eq_path, md_prod_path)
        pdb_files, itp_files, top_files = process_all_topologies(m_smiles, mols, a_smiles, ans, c_smiles, cats, ions_pdb, ff_path, pack_path, md_em_path, md_eq_path, md_prod_path)
        print("\nDone!\n")

        print("#################################")
        print("\nCreating electrolyte structure ...")
        print("\n#################################\n")  

        a_side, n_mols_box = number_of_molecules(m_conc, a_conc, c_conc, i_conc)
        system, packmol_file = packmol_build(work_path, pack_path, md_em_path, pdb_files, a_side, n_mols_box, ions)
        
        print(f"\nDone! Check now the {packmol_file}! \n")

        print("#################################")
        print("\nPreparing  simulations ...")
        print("\n#################################\n")  

        print("Adjusting topologies...")
        
        prepare_topol(doc_id, mols, ans, cats, ions, n_mols_box, md_em_path, md_eq_path, md_prod_path, config)
        
        print("\nDone!\n")
        
        print("#################################")
        print("\nRunning MD simulations ...")
        print("\n#################################\n") 
        
        md_simulation_run(BASE_DIR, config, work_path, md_em_path, md_eq_path, md_prod_path)
        
        print("\nDone!\n")
        
        

if __name__ == "__main__":
    main()