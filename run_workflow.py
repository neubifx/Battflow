import os

from battflow.database import db_connection, scan_properties_collection
from battflow.md_setup import smiles_setup, prepare_simulation_paths, prepare_molecule_topologies, prepare_anion_topologies

def main():
    print("Connecting to DB ...")
    db, collection = db_connection()

    print("Scanning collections for missing properties ...")
    flag, doc_id = scan_properties_collection(collection)
    
    if flag:
        print(f"Found missing property in document ID {doc_id}!")
        doc = collection.find_one({"_id" : doc_id})
        print("Building up electrolyte structure ...")
        mols, ans, m_smiles, a_smiles = smiles_setup(doc)
        print("Molecules:", mols)
        print("Anions:", ans)
        print("Setting up simulation folders...")
        work_path, setup_path = prepare_simulation_paths(doc_id)
        print("Done!\n")
        print("#################################")
        print("\nPreparing molecules topologies ...")
        print("\n#################################\n")
        prepare_molecule_topologies(work_path, setup_path, mols, m_smiles)
        print("\nDone!\n")
        print("#################################")
        print("\nPreparing anions topologies ...")
        print("\n#################################\n")
        prepare_anion_topologies(work_path, setup_path, ans, a_smiles, m_smiles)
        print("\nDone!\n")        
        

if __name__ == "__main__":
    main()