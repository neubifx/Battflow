from battflow.database import db_connection, scan_properties_collection
from battflow.md_setup import smiles_setup

def main():
    print("Connecting to DB ...")
    db, collection = db_connection()

    print("Scanning collections for missing properties ...")
    flag, doc_id = scan_properties_collection(collection = collection)
    
    if flag:
        print(f"Found missing property in document ID {doc_id}!")
        doc = collection.find_one({"_id" : doc_id})
        print("Building up electrolyte structure ...")
        mols_smiles, anions_smiles = smiles_setup(doc = doc)
        print("Molecules:", mols_smiles)
        print("Anions:", anions_smiles)

if __name__ == "__main__":
    main()