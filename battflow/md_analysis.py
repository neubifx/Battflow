from pathlib import Path
import os
import re
import numpy as np
from scipy.stats import linregress
#from battflow.database import db_connection

import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
from solvation_analysis.solute import Solute

from bson import ObjectId

from battflow.utils import normalize_id

def setup_mda_analysis(md_prod_path, mols, ans, solute_ion = "li"):
    """
    Set up the MDAnalysis universe and atom selections for trajectory analysis.

    Parameters
    ----------
    md_prod_path : pathlib.Path
        Directory containing the production simulation topology (.tpr)
        and trajectory (.xtc) files.
    mols : list of str
        Names of the molecular components in the electrolyte.
    ans : list of str
        Names of the anionic components in the electrolyte.
    solute_ion : str, optional
        Elemental symbol of the solute ion. Supported ions are Li, Na,
        K, Ca and Zn. Default is "li".

    Returns
    -------
    u : MDAnalysis.Universe
        Universe containing the production simulation trajectory.
    mda_names : list of str
        Names of the molecular and anionic components.
    mda_resnames : set of str
        Unique residue names present in the simulation system.
    dict_solvation : dict
        Mapping of component names to MDAnalysis AtomGroup objects.
    ion_solute : MDAnalysis.AtomGroup
        Atom selection corresponding to the solute ion.

    Raises
    ------
    ValueError
        If the specified solute ion is not supported.

    Notes
    -----
    Molecular and anionic components are selected using their sequential
    numerical residue names. Additional cationic components are excluded
    from the solvation analysis.
    """
    
    u = mda.Universe(md_prod_path / "prod_0_1.tpr", md_prod_path / "prod_0_1.xtc") 
    mda_resnames = set(u.residues.resnames)
    mda_names = mols + ans # cations are not included in the solvation structure analysis and transference number calcs

    #variables setup
    dict_solvation = {}
    var_names = []
    for i, name in enumerate(mda_names, start =1):
        var_name = f"c_{i:03d}"
        var_names.append(var_name)
                #use flexible variable depending on the number of components of the electrolyte
        globals()[var_name] = u.select_atoms(f"resname {i:03d}")
        
        #create dict to use in Solvation Analysis
        dict_solvation[name] = globals()[var_name]

    #Flexible definition of the solute
    solute_ion = solute_ion.lower()

    #create dict to change resnames accordingly to the args

    ion_resnames = {
        "li": "LI",
        "k":  "K",
        "na": "NA",
        "ca": "CA",
        "zn": "ZN",
    }

    if solute_ion not in ion_resnames:
        raise ValueError(f"Unsupported solute ion: {solute_ion}")

    ion_resname = ion_resnames[solute_ion]

    #select atoms usiing mduniverse
    ion_solute = u.select_atoms(f"resname {ion_resname}")
        
    return u, mda_names, mda_resnames, dict_solvation, ion_solute
    
def solvation_structure_analysis(u, ion_solute, dict_solvation, solute_ion = "li"):
    """
    Calculate the solvation structure of the selected solute ion.

    The analysis is performed over the second half of the production
    trajectory using the solvation_analysis package.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Universe containing the production simulation trajectory.
    ion_solute : MDAnalysis.AtomGroup
        Atom selection corresponding to the solute ion.
    dict_solvation : dict
        Mapping of electrolyte component names to MDAnalysis AtomGroup
        objects.
    solute_ion : str, optional
        Elemental symbol of the solute ion. Default is "li".

    Returns
    -------
    solute : solvation_analysis.solute.Solute
        Solute object containing the completed solvation analysis.
    coordination_number : dict
        Average coordination numbers of the electrolyte components
        around the solute ion.
    pairing_percentage : dict
        Pairing fractions between the solute ion and the electrolyte
        components.
    solvation_shell : list of dict
        Relative populations of the three most common solvation
        structures.

    Notes
    -----
    Solvation radii are initially determined automatically. If a radius
    cannot be identified for a component and its name is recognised in
    the error message, a default radius of 3.0 angstroms is assigned
    to that component and the analysis is repeated.
    """
    n_frames = len(u.trajectory)

    # use solute_ion in uppercase name the solute
    solute_name = solute_ion.upper()

    try:
        # First attempt without specifying radii
        solute = Solute.from_atoms(ion_solute, dict_solvation, solute_name = solute_name)
        solute.run(start=int(n_frames / 2), stop=n_frames, step=1)

    except AssertionError as e:
        msg = str(e)
        missing_radii = {}

        if "could not identify a solvation radius for" in msg:
            # Extract everything after 'for', remove punctuation
            after_for = msg.split("for", 1)[-1]
            clean_text = after_for.replace(".", " ").replace(",", " ")
            # Split into words and filter by dict_solvation keys
            words = clean_text.split()
            for word in words:
                if word in dict_solvation:
                    missing_radii[word] = 3.0
        else:
            raise  # If it's a different AssertionError, re-raise it

        if missing_radii:
            print(f"[INFO] Assigning default radius=3.0 for solvents: {', '.join(missing_radii.keys())}")
        else:
            print("[WARN] No solvents matched from error message — radius defaults not applied.")

        # Retry with missing radii
        solute = Solute.from_atoms(
            ion_solute,
            dict_solvation,
            solute_name= solute_name,
            radii=missing_radii
        )
        solute.run(start=int(n_frames / 2), stop=n_frames, step=1)

    # Collect results
    coordination_number = solute.coordination.coordination_numbers
    pairing_percentage = solute.pairing.solvent_pairing
    solvation_shell = solute.speciation.speciation_fraction.head(3).to_dict(orient="records")

    return solute, coordination_number, pairing_percentage, solvation_shell

def get_diffusion(u, selection):
    """
    Calculate the self-diffusion coefficient of a selected group of atoms.

    The diffusion coefficient is obtained from the slope of the
    three-dimensional mean squared displacement (MSD) using the
    Einstein relation.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Universe containing the production simulation trajectory.
    selection : str
        MDAnalysis atom selection string identifying the atoms for
        which the diffusion coefficient is calculated.

    Returns
    -------
    D : float
        Self-diffusion coefficient in angstroms squared per picosecond
        (angstrom^2/ps).

    Notes
    -----
    The MSD is calculated using the FFT-based EinsteinMSD algorithm.
    Linear regression is performed over the entire available MSD
    time series. The diffusion coefficient is calculated as the
    fitted slope divided by six.

    The trajectory must contain unwrapped coordinates to avoid
    discontinuities arising from periodic boundary conditions.
    """
    # Compute MSD
    msd_calc = EinsteinMSD(u, select=selection, msd_type="xyz", fft=True)
    msd_calc.run()

    msd_arr = msd_calc.results.timeseries     # Å²
    dt = u.trajectory.dt                      # ps
    t = np.arange(len(msd_arr)) * dt          # ps

    n_frames = len(u.trajectory)
    i0 = 0
    i1 = n_frames

    slope = linregress(t[i0:i1], msd_arr[i0:i1]).slope
    D = slope / (2 * 3)                        # Å²/ps | 2 * 3(xyz)
    D_m2_s = D * 1e-20 / 1e-12                 # m²/s not using right now

    return D

def ions_anions_transference_number(u, mols, ans, a_conc, ions, i_conc, solute_ion = "li"):
    """
    Calculate ionic diffusion coefficients and estimate the solute-ion
    transference number.

    The transference number is estimated from the self-diffusion
    coefficients and the number of ions of each species.

    Parameters
    ----------
    u : MDAnalysis.Universe
        Universe containing the production simulation trajectory.
    mols : list of str
        Names of the molecular components in the electrolyte.
    ans : list of str
        Names of the anionic components.
    a_conc : list of str
        Molar concentrations of the anionic components, expressed as
        strings with the suffix "M".
    ions : list of str
        Elemental symbols of the solute ions.
    i_conc : list of str
        Molar concentrations of the solute ions, expressed as strings
        with the suffix "M".
    solute_ion : str, optional
        Elemental symbol of the solute ion. Default is "li".

    Returns
    -------
    D_solute : float
        Self-diffusion coefficient of the solute ion in angstrom^2/ps.
    D_ans_dict : dict
        Mapping of anion names to their self-diffusion coefficients
        in angstrom^2/ps.
    t_final : float
        Estimated dimensionless solute-ion transference number.

    Raises
    ------
    ValueError
        If the specified solute ion is not supported.

    Notes
    -----
    The number of ions is estimated from the input concentrations
    using a cubic box with a side length of 50 angstroms.

    The current implementation determines the solute-ion population
    from the lithium concentration. Additional cationic components,
    ionic charge weighting and correlated ionic motion are not
    included in the transference number calculation.
    """

    solute_ion = solute_ion.lower()

    ion_resnames = {
        "li": "LI",
        "k":  "K",
        "na": "NA",
        "ca": "CA",
        "zn": "ZN",
    }

    if solute_ion not in ion_resnames:
        raise ValueError(f"Unsupported solute ion: {solute_ion}")

    solute_resname = ion_resnames[solute_ion]

    #box settings 
    a_side = 50 # fixed 50 Angstroms box
    box_vol = (a_side**3) * 1e-27 # box volume in L
    avo = 6.022e23
             
    #calculate number of solute ions
    packed_ions_concs = list(zip(ions, i_conc))
    
    for ion, conc in packed_ions_concs:
        if ion == "li":
            conc_value = float(conc.replace("M", ""))
            i_mol = round(conc_value * box_vol * avo)

    #Diff coeff and transf number for the solute
    D_solute = get_diffusion(u, f"resname {solute_resname}")
    t_n_solute = D_solute*i_mol

    #calculate number of anions
    a_mols_box = []
    for conc in a_conc:
        conc_value = float(conc.replace("M", "")) #removes concentration from the string and turn it into a value
        a_mol = round(conc_value * box_vol * avo)
        a_mols_box.append(a_mol)  
    
    packed_ans_concs = list(zip(ans, a_mols_box))

    #Transference number and diffusion coefficients
    t_n_ans = []
    D_ans = []
    start = len(mols)+1
    for i, (a, a_mol) in enumerate(packed_ans_concs, start = start):
        D_an = get_diffusion(u, f"resname {i:03d}") #i starts at the end of the mols
        D_ans.append(D_an)
        t_n_an = D_an*a_mol
        t_n_ans.append(t_n_an)

    #calculate transference number
    t_n_ans.append(t_n_solute)
    t_final = t_n_solute / sum(t_n_ans)

    # Make dict for ans
    D_ans_dict = {a: float(D) for a, D in zip(ans, D_ans)}

    return D_solute, D_ans_dict, t_final

from bson import ObjectId

def upload_calculated_data(collection, doc_id, 
                           coordination_number, pairing_percentage, 
                           solvation_shell, D_solute, D_ans_dict, t_final):
    """
    Upload calculated molecular dynamics properties to MongoDB.

    Parameters
    ----------
    collection : pymongo.collection.Collection
        MongoDB collection containing the electrolyte documents.
    doc_id : str or bson.ObjectId
        Identifier of the electrolyte document to be updated.
    coordination_number : dict
        Calculated coordination numbers of the electrolyte components.
    pairing_percentage : dict
        Calculated pairing fractions for the solute ion.
    solvation_shell : list of dict
        Statistics of the most common solvation structures.
    D_solute : float
        Self-diffusion coefficient of the solute ion in angstrom^2/ps.
    D_ans_dict : dict
        Anion names and their corresponding self-diffusion coefficients
        in angstrom^2/ps.
    t_final : float
        Estimated dimensionless solute-ion transference number.

    Returns
    -------
    None

    Notes
    -----
    The document identifier is normalised before the database update.
    The calculated properties are stored within the simulation_data
    field of the corresponding electrolyte document.
    """
    
    safe_id = normalize_id(doc_id)

    collection.update_one(
        {"_id": safe_id},
        {"$set": {
            "simulation_data.coordination_number": coordination_number,
            "simulation_data.pairing_percentage": pairing_percentage,
            "simulation_data.solvation_statistics": solvation_shell,
            "simulation_data.diffusion_coefficients": {
                "ions": D_solute,
                "anions": D_ans_dict
            },
            "simulation_data.transference_number": t_final
        }}
    )

