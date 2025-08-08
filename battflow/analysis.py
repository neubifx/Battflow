from pathlib import Path
import os

import MDAnalysis as mda
from MDAnalysis.analysis.msd import EinsteinMSD
from solvation_analysis.solute import Solute

from bson import ObjectId

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
    
    solute = Solute.from_atoms(ion_solute, dict_solvation, solute_name="Li")
    #running analysis
    n_frames = len(u.trajectory)
    solute.run(start = int(n_frames/2), stop = n_frames, step = 1)

    # inspect the coordination numbers. Interpreted as the mean number of solvent coordinated with each solute
    coordination_number = solute.coordination.coordination_numbers

    # inspect the pairing percentages. The percentage of solutes paired with each solvent
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

def upload_calculated_data(doc_id, coordination_number, pairing_percentage, solvation_shell):
    # Update both fields in one go
    collection.update_one(
        {"_id": ObjectId(doc_id)},  # document to update
        {"$set": {
            "properties.coordination_number": coordination_number,
            "properties.pairing_percentage": pairing_percentage,
            "properties.solvation_statistics": solvation_shell,
            "properties.diffusion_coefficients.ions": D_solute,
            "properties.diffusion_coefficients.anions": D_ans_dict,
            "properties.transference_number": t_final
        }}
    )