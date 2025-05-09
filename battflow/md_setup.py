
def smiles_setup(doc):
    """
    Retrieves and processes molecular components and their SMILES
    representations from a MongoDB database.
    
    Returns:
        A list of tuples containing (component, SMILES) pairs or None if no 
        document is found.
    """
            
    #electrolyte composition
    molecules = doc["components"]["molecules"]
    anions = doc["components"]["anions"]
    
    #smiles
    smiles_molecules = doc["smiles"]["molecules"]
    smiles_anions = doc["smiles"]["anions"]

    mols_smiles = list(zip(molecules, smiles_molecules))
    anions_smiles = list(zip(anions, smiles_anions))
    
    return mols_smiles, anions_smiles
