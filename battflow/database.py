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
    Scan over the the documents in the collection to find if there are properties
    to be calculated. The function will stop in the first missing property found and
    the _id of the document will be stored for further calculation.

    Returns:
        tuple:
            - bool: True if a missing property is found. False otherwise.
            - ObjectID or None: The _id of the document with missing properties as a MongoDB's ObjectID,
              or None if no missing property is found.
    """
    
    for doc in collection.find():
        print(f"Checking document {doc['_id']} ...")
        
        #add a flag if found
        if dict_null_check(doc["properties"]):
            doc_id = doc["_id"]
            return True, doc_id # returning both. Remember to unpack tuple later
            
    print("There are no properties to calculate at this moment")    
    return False, None #getting out of the if loop to be reachable



