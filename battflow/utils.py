from bson import ObjectId

def normalize_id(doc_id):
    """
    Safely normalize the document ID.
    - If doc_id is already an ObjectId, return it.
    - If doc_id is a string that looks like a valid ObjectId, convert it.
    - Otherwise, treat the ID as a plain string (DOI, custom key, etc.).
    """
    if isinstance(doc_id, ObjectId):
        return doc_id

    if isinstance(doc_id, str) and ObjectId.is_valid(doc_id):
        return ObjectId(doc_id)

    # Otherwise, it is something like a DOI / custom string ID
    return doc_id