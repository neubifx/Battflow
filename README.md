# Battflow - Automated Workflow for Properties Prediction of Liquid Electrolytes

Battflow is an automated workflow designed to identify and compute missing molecular dynamics (MD) and density functional theory (DFT) properties in a MongoDB-based liquid electrolyte database.

## Overview

The workflow is composed of the following steps:
1. **Query MongoDB** for entries with missing properties.
2. **Preprocess** structures from SMILES information and setup force-fields and electrolyte composition.
3. **Run MD simulations** using GROMACS to compute diffusion properties and solvation-cluster statistics.
4. **Perform molecular calculations** using ORCA to obtain binding energies, HOMO, LUMO, and HOMO–LUMO gaps for individual components and the most populated solvation clusters.
5. **Postprocess** results and update MongoDB.

## Installation

1. Create the Conda environment

This workflow relies heavily on ACPYPE (https://github.com/alanwilter/acpype) for GAFF2 topology generation. Therefore, an Anaconda environment and a Linux-based system are required for successful execution.

```
conda env create -f environment.yml
conda activate battflow_env
```

2. External software requirements

The workflow requires the following external software:

- Python ≥ 3.9
- [MongoDB](https://www.mongodb.com/)
- [ASE](https://wiki.fysik.dtu.dk/ase/)
- [GROMACS](https://www.gromacs.org/)
- [ORCA](https://orcaforum.kofo.mpg.de/)

3. Configuration

Ensure that `gmx`, `orca`, and any MPI executables are available in your system `PATH` when the Battflow environment is active. All configuration parameters are defined in:

```
config/default.yaml
```

This file can be modified to adapt the workflow to local machines or HPC environments.

### MongoDB configuration

```
mongodb:
    host: localhost                   # Replace with your MongoDB Atlas host if preferred
    port: 27017
    database: working_db
    collection: working_collection
```

### MD run environment

```
md_run_env:
    mdrun: "gmx mdrun"     # Command to run GROMACS
    mpiexec: ""            # Path to the MPI executor (e.g. output of "which srun")
    ncores: 8
```

The directory `resources/md_run_files/` contains default GROMACS input files for minimisation, equilibration, and production runs. 

### DFT simulation configuration

```
dft_simulations:
    orca_profile: /usr/local/orca_5.0.3/orca    # Full path to the orca executable
    orca_input_block: "B3LYP 6-311+G(d,p)"      # You can change the level of theory and add additional parameters for ORCA simulations
    ncores: 8                                   # Number of cores
    li_energy: -203.567718                      # Energy of Li atom calculated at the B3LYP 6-311+G(d,p) level of theory. This value should match the methodology adopted on your simulations
```

## Running the workflow

Typical usage:

```
python run_workflow.py --config config/default.yaml
```

The workflow extracts SMILES and concentration information from documents stored in MongoDB.
It expects documents with the following structure:


```
{
  "_id": { "$oid": "681ba5b9aabcaf427f287f62" },
  "components": {
    "molecules": ["ec", "emc"],
    "anions": ["pf6"],
    "cations": ["methylimidazolium"],
    "ions": ["li"]
  },
  "smiles": {
    "molecules": ["C1COC(=O)O1", "CCOC(=O)C"],
    "anions": ["F[P-](F)(F)(F)(F)F"],
    "cations": ["C[N+]1=CNC=C1"]
  },
  "concentrations": {
    "molecules": ["1M", "0.8M"],
    "anions": ["1.2M"],
    "cations": ["0.2"],
    "ions": ["1M"]
  },
  "simulation_data": {
    "diffusion_coefficients": { "ions": null, "anions": null },
    "transference_number": null,
    "solvation_statistics": null,
    "coordination_number": null,
    "pairing_percentage": null,
    "dft_energies": null
  }
}


```

### Output fields populated by Battflow

1. Molecular dynamics transport and solvation properties

The field `simulation_data` will be filled with:

- Diffusion coefficients,
- Transference number,
- The three most common solvation structures,
- Their Cartesian coordinates (XYZ format),
- Coordination number for the solvated ion,
- Pairing percentage for the solvated ion


2. DFT energies of individual components

The field `dft_energies` will store:

- DFT energy,
- HOMO,
- LUMO,
- HOMO–LUMO gap

for each molecular component in the electrolyte formulation.








