import os
import argparse
import logging
import traceback
from pathlib import Path
import yaml

from battflow.database import db_connection, scan_properties_collection, config_path
from battflow.md_setup import (
    names_smiles_molarity_setup,
    prepare_simulation_paths,
    prepare_molecule_topologies,
    prepare_anion_topologies,
    prepare_cation_topologies,
    process_ion_topologies,
    process_all_topologies,
    number_of_molecules,
    packmol_build,
)
from battflow.topol_tools import prepare_topol
from battflow.md_run import md_simulation_run

from battflow.topol_tools import prepare_topol
from battflow.md_run import md_simulation_run
from battflow.md_analysis import (
    setup_mda_analysis,
    solvation_structure_analysis,
    get_diffusion,
    ions_anions_transference_number,
    upload_calculated_data
)

from battflow.dft_analysis import (
    create_dft_folders_and_coordinates,
    extract_solvation_structure_from_md,
    extract_orbital_energies,
    calculate_component_dft_energies,
    upload_dft_calculated_data   
)


def main():
    parser = argparse.ArgumentParser(description="Run Battflow workflow.")
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to alternative YAML config file (default: config/default.yaml)",
    )
    args = parser.parse_args()

    # Load YAML config (defaults to package's default if none provided)
    BASE_DIR, config = config_path(args.config)

    # Set working base directory to current working directory
    WORK_DIR = Path.cwd()

    print("Connecting to DB ...")
    db, collection = db_connection(config)

    print("Scanning collections for missing properties ...")
    flag, doc_ids = scan_properties_collection(collection)

    if flag:
        failed = []
        succeeded = []

        for doc_id in doc_ids:
            print(f"Found missing property in document ID {doc_id}!")
            try:
                doc = collection.find_one({"_id": doc_id})

                print("Building up electrolyte structure ...")
                mols, ans, cats, ions, m_smiles, a_smiles, c_smiles, m_conc, a_conc, c_conc, i_conc = names_smiles_molarity_setup(doc)

                print("Molecules:", mols)
                print("Anions:", ans)
                print("Cations:", cats)
                print("Ions:", ions)
                print("Setting up simulation folders...")

                # Use WORK_DIR for all folder creation
                paths = prepare_simulation_paths(doc_id, base_dir=WORK_DIR)
                (work_path, setup_path, ff_path, pack_path, md_path,
                 md_em_path, md_eq_path, md_prod_path, dft_path) = paths

                print("Done!\n")
                print("#################################")
                print("\nPreparing molecules topologies ...")
                print("\n#################################\n")

                prepare_molecule_topologies(work_path, ff_path, mols, m_smiles)

                print("\nDone!\n")
                print("#################################")
                print("\nPreparing anions topologies ...")
                print("\n#################################\n")

                prepare_anion_topologies(work_path, ff_path, ans, a_smiles, m_smiles)

                print("\nDone!\n")
                print("#################################")
                print("\nPreparing cations topologies ...")
                print("\n#################################\n")

                prepare_cation_topologies(work_path, ff_path, cats, c_smiles, a_smiles, m_smiles)

                print("\nDone!\n")

                print("Processing topologies files ...")
                ions_itp_file, topol_main_file, ions_pdb = process_ion_topologies(
                    BASE_DIR, config, ions, pack_path, md_em_path, md_eq_path, md_prod_path)
                pdb_files, itp_files, top_files = process_all_topologies(
                    m_smiles, mols, a_smiles, ans, c_smiles, cats, ions_pdb,
                    ff_path, pack_path, md_em_path, md_eq_path, md_prod_path)
                print("\nDone!\n")

                print("#################################")
                print("\nCreating electrolyte structure ...")
                print("\n#################################\n")

                a_side, n_mols_box = number_of_molecules(m_conc, a_conc, c_conc, i_conc)
                system, packmol_file = packmol_build(work_path, pack_path, md_em_path, pdb_files, a_side, n_mols_box, ions)

                print(f"\nDone! Check now the {packmol_file}! \n")

                print("#################################")
                print("\nPreparing simulations ...")
                print("\n#################################\n")

                print("Adjusting topologies...")
                prepare_topol(doc_id, mols, ans, cats, ions, n_mols_box, md_em_path, md_eq_path, md_prod_path, config)

                print("\nDone!\n")

                print("#################################")
                print("\nRunning MD simulations ...")
                print("\n#################################\n")

                md_simulation_run(BASE_DIR, config, work_path, md_em_path, md_eq_path, md_prod_path)

                print("\nDone!\n")
                succeeded.append(doc_id)
                logging.info("Document %s finished successfully", doc_id)
                print(f"\nDocument {doc_id} completed successfully!\n")

                print("#################################")
                print("\nRunning Analysis from MD simulations ...")
                print("\n#################################\n") 
                
                u, mda_names, mda_resnames, dict_solvation, ion_solute = setup_mda_analysis(md_prod_path, mols, ans)
                
                solute, coordination_number, pairing_percentage, solvation_shell = solvation_structure_analysis(u, ion_solute, dict_solvation)
                
                D_solute, D_ans_dict, t_final = ions_anions_transference_number(u, mols, ans, a_conc, ions, i_conc)
                
                print("#################################")
                print("\nUploading MD data into MongoDB ...")
                print("\n#################################\n") 
                               
                upload_calculated_data(collection, doc_id, coordination_number, pairing_percentage, solvation_shell, D_solute, D_ans_dict, t_final)
                
                print("\nDone!\n")
                
                print("#################################")
                print("\nSetting up DFT simulations ...")
                print("\n#################################\n")
                
                
                dft_component_folders, packed_molecules, packed_anions, packed_cations = create_dft_folders_and_coordinates(m_smiles, mols, a_smiles, ans, c_smiles, cats, ions_pdb, dft_path, pdb_files)
                
                print("\nDone!\n")
                
                print("\nExtracting solvation structures from MD simulations\n")
                
                solvation_shell = extract_solvation_structure_from_md(u, solvation_shell, solute)
                
                print("\nDone!\n")
                
                print("#################################")
                print("\n Running DFT simulations ...")
                print("\n#################################\n")  

                print("\nStarting DFT calculations for each electrolyte component\n")                
 
                energies_dict = calculate_component_dft_energies(dft_component_folders, mols, ans, cats, config, packed_molecules, packed_cations, packed_anions)
                
                print("\nDone!\n")
                
                print("#################################")
                print("\nUploading MD data into MongoDB ...")
                print("\n#################################\n") 
                
                upload_dft_calculated_data(doc_id, solvation_shell, energies_dict, collection)
                
                print("\nDone!\n")                
                

            except Exception as e:
                failed.append(doc_id)
                print(f"\n>>> ERROR with document {doc_id}: {e}\n")
                logging.error("Document %s failed with error: %s", doc_id, e)
                logging.error(traceback.format_exc())

        print("#################################")
        print("Summary of runs:")
        print(f"  Successes: {len(succeeded)}")
        print(f"  Failures : {len(failed)}")
        if failed:
            print("  Failed IDs:", failed)
        print("#################################")

        logging.info("Summary: %d successes, %d failures", len(succeeded), len(failed))
        if failed:
            logging.info("Failed IDs: %s", failed)


if __name__ == "__main__":
    # configure logging at the start
    logging.basicConfig(
        filename="md_run.log",
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    main()

