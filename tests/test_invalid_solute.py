# tests/test_invalid_solute.py

import pytest

def validate_solute(solute_ion):

    ion_resnames = {
        "li": "LI",
        "k": "K",
        "na": "NA",
        "ca": "CA",
        "zn": "ZN",
    }

    solute_ion = solute_ion.lower()

    if solute_ion not in ion_resnames:
        raise ValueError(f"Unsupported solute ion: {solute_ion}")

def test_invalid_solute():

    with pytest.raises(ValueError):
        validate_solute("mg")
