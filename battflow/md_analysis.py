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

from Battflow.battflow.utils import normalize_id

def setup_mda_analysis(md_prod_path, mols, ans):
    
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

    #rigid definition of Li as a solute. Possibly changing to others in the submission arguments
    ion_solute = u.select_atoms("resname LI")    
        
    return u, mda_names, mda_resnames, dict_solvation, ion_solute
    
def solvation_structure_analysis(u, ion_solute, dict_solvation):
    n_frames = len(u.trajectory)

    try:
        # First attempt without specifying radii
        solute = Solute.from_atoms(ion_solute, dict_solvation, solute_name="Li")
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
            solute_name="Li",
            radii=missing_radii
        )
        solute.run(start=int(n_frames / 2), stop=n_frames, step=1)

    # Collect results
    coordination_number = solute.coordination.coordination_numbers
    pairing_percentage = solute.pairing.solvent_pairing
    solvation_shell = solute.speciation.speciation_fraction.head(3).to_dict(orient="records")

    return solute, coordination_number, pairing_percentage, solvation_shell

def get_diffusion(u, selection):
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

def ions_anions_transference_number(u, mols, ans, a_conc, ions, i_conc):

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
    D_solute = get_diffusion(u, "resname LI")
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
    
    safe_id = normalize_id(doc_id)

    collection.update_one(
        {"_id": safe_id},
        [
            {
                "$set": {
                    "simulation_data.diffusion_coefficients": {
                        "$ifNull": ["$simulation_data.diffusion_coefficients", {}]
                    }
                }
            },
            {
                "$set": {
                    "simulation_data.coordination_number": coordination_number,
                    "simulation_data.pairing_percentage": pairing_percentage,
                    "simulation_data.solvation_statistics": solvation_shell,
                    "simulation_data.diffusion_coefficients.ions": D_solute,
                    "simulation_data.diffusion_coefficients.anions": D_ans_dict,
                    "simulation_data.transference_number": t_final
                }
            }
        ]
    )
