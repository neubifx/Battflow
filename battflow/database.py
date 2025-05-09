from pymongo import MongoClient



def db_connection():
    """
    Temporary function to connect with a local deployment of MongoDB
    """
    client = MongoClient("mongodb://localhost:27017/")

    db = client["working_db"]
    collection = db["working_collection"]
    return db, collection
    
def dict_null_check(data, path = ""):
    """
    Recursively scan a nested dictionary in a single MongoDB document for None values.

    Parameters:
        data (dict):  
            The dictionary path to scan.  

        path (str):  
            The “dotted” key path of the current recursion level, to printout 
            the location in the nested dict . 

    Returns:
        bool:  
            True if any value in the dictionary (or its nested dicts) is None,
            False otherwise. Used to flag subsequent simulations.
    """
    
    for key, value in data.items():
        parent_path = f"{path}.{key}" if path !="" else key 
        
        if value is None:
            print(f"{parent_path} is None")
            return True
            
        elif isinstance(value, dict): 
            if dict_null_check(value, path = parent_path):
                return True
                
    return False

def scan_properties_collection():
    """
    Scan over the the documents in the collection to find if there are properties
    to be calculated. The function will stop in the first missing property found and
    the _id of the document will be stored for further calculation.

    Returns:
        tuple:
            - bool: True if a missing property is found. False otherwise.
            - ObjectID or None: The _id of the document with missing properties as a MongoDB's ObjectID,
              or None if no missing property is found.
    """
    db, collection = db_connection()
    
    for doc in collection.find():
        print(f"Checking document {doc['_id']} ...")
        
        #add a flag if found
        if dict_null_check(doc["properties"]):
            print(f"Found properties to calculate in {doc['_id']}")
            doc_id = doc["_id"]
            return True, doc_id # returning both. Remember to unpack tuple later
            
        else:
            print("There are no properties to calculate at this moment")
            return False, None

def smiles_setup():
    """
    Retrieves and processes molecular components and their SMILES
    representations from a MongoDB database. It is based on the document
    output from scan_properties_collection()
    
    Returns:
        A list of tuples containing (component, SMILES) pairs or None if no 
        document is found.
    """
    flag, doc_id = scan_properties_collection();

    if flag:
        doc = collection.find_one({"_id": doc_id})
            
        #electrolyte composition
        molecules = doc["components"]["molecules"]
        anions = doc["components"]["anions"]
        
        #smiles
        smiles_molecules = doc["smiles"]["molecules"]
        smiles_anions = doc["smiles"]["anions"]
    
        mols_smiles = list(zip(molecules, smiles_molecules))
        anions_smiles = list(zip(anions, smiles_anions))
    return mols_smiles, anions_smiles

