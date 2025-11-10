import yaml
from pathlib import Path

from pymongo import MongoClient

def config_path(config_file_path=None):
    """
    Point towards the .yaml file containing diverse settings

    Args:
        config_file_path (str or Path, optional): Path to the config YAML file. If None, use default.

    Returns:
       BASE_DIR (pathlib.Path): Base directory of Battflow
       config (dict): Dictionary containing the loaded .yaml file 
    """
    BASE_DIR = Path(__file__).resolve().parents[1]
    if config_file_path is not None:
        config_file = Path(config_file_path)
    else:
        config_file = BASE_DIR / "config" / "default.yaml"
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    return BASE_DIR, config


def db_connection(config):
    """
    Function to connect with a local deployment of MongoDB
    """
    client = MongoClient(config["mongodb"]["host"])

    db = client[config["mongodb"]["database"]]
    collection = db[config["mongodb"]["collection"]]
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
        parent_path = f"{path}.{key}" if path !="" else key # Can return parent_path if spitting out proeprty is needed
        
        if value is None:
            return True
            
        elif isinstance(value, dict): 
            if dict_null_check(value, path = parent_path):
                return True
                
    return False

def scan_properties_collection(collection):
    """
    Scan over the documents in the collection to find if there are properties
    to be calculated. Returns a list of _id for all documents with missing properties.

    Returns:
        tuple:
            - bool: True if any missing property is found. False otherwise.
            - list: List of _id for documents with missing properties.
    """
    doc_ids = []
    for doc in collection.find():
        print(f"Checking document {doc['_id']} ...")
        
        if dict_null_check(doc["simulation_data"]):
            doc_ids.append(doc["_id"])
            
    if doc_ids:
        print(f"\n{len(doc_ids)} documents with missing properties found.")
        return True, doc_ids
        
    print("There are no properties to calculate at this moment")
    return False, []



