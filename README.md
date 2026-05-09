![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Version](https://img.shields.io/badge/version-1.0.0-orange)

# Battflow - Automated Workflow for Property Prediction of Liquid Electrolytes

Battflow is an automated workflow designed to identify and compute battery electrolyte properties using molecular dynamics (MD) and density functional theory (DFT).

## Overview

The workflow is composed of the following steps:

1. **Query MongoDB** for entries with missing properties.
2. **Preprocess** structures from SMILES information and set up force fields and electrolyte compositions.
3. **Run MD simulations** using GROMACS to compute diffusion properties and solvation-cluster statistics.
4. **Perform molecular calculations** using ORCA to obtain binding energies, HOMO, LUMO, and HOMO–LUMO gaps for individual components and the most populated solvation clusters.
5. **Postprocess** the results and update MongoDB.

---

## Installation

### 1. Install Battflow

Battflow relies on [ACPYPE](https://github.com/alanwilter/acpype) for GAFF2 topology generation. Therefore, a Linux-based system is recommended for successful installation and execution.

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate battflow_env
```

Install Battflow:

```bash
pip install git+https://github.com/neubifx/Battflow.git
```

For editable/development installation:

```bash
git clone https://github.com/neubifx/Battflow.git
cd Battflow
pip install -e .
```

---

### 2. External software requirements

Battflow interfaces with external software packages that must be available in the system environment:

- Python ≥ 3.10
- [MongoDB](https://www.mongodb.com/)
- [GROMACS](https://www.gromacs.org/)
- [ORCA](https://orcaforum.kofo.mpg.de/)

In particular, `gmx`, `orca`, and any required MPI executables should be accessible through the system `PATH` when the Battflow environment is active.

---

## Configuration

Users can inspect and modify the default configuration according to their simulation environment. Generate a copy of the default configuration file in the current working directory:

```bash
battflow --write-config
```

This creates:

```text
default.yaml
```

Users can then modify the configuration file according to their local environment (e.g., MongoDB URI, ORCA and GROMACS commands, HPC submission commands, number of cores, etc.).

---

### MongoDB configuration

This section should point to the database and collection from which Battflow will retrieve electrolyte compositions for calculation. The default options correspond to a local MongoDB installation, although a MongoDB Atlas URI can also be provided in `host`.

```yaml
mongodb:
    host: localhost                   # Replace with your MongoDB Atlas host if preferred
    port: 27017
    database: working_db
    collection: working_collection
```

---

### MD run environment

The default options assume a standard Linux installation of GROMACS. These commands can be adapted to match the commands used in an HPC environment. The number of cores used during workflow execution can also be adjusted.

```yaml
md_run_env:
    mdrun: "gmx mdrun"     # Command used to run GROMACS
    mpiexec: ""            # Path to the MPI executor (e.g. output of "which srun")
    ncores: 8
```

---

### DFT simulation configuration

The default options assume a standard Linux installation of ORCA. These commands can be adapted to match the commands used in an HPC environment. The number of cores used during workflow execution can also be adjusted.

```yaml
dft_simulations:
    orca_profile: /usr/local/orca_5.0.3/orca    # Full path to the ORCA executable
    orca_input_block: "B3LYP 6-311+G(d,p)"      # Level of theory and additional ORCA parameters
    ncores: 8                                   # Number of cores
    li_energy: -203.567718                      # Energy of Li atom calculated at the selected level of theory
```

---

## Running Battflow

Run Battflow using the configuration file:

```bash
battflow --config default.yaml
```

Optional: specify the solute ion used during MD analysis. If no option is provided, Li will be considered the solute ion.

```bash
battflow --config default.yaml --solute-ion li
```

Supported solute ions currently include:

- Li
- Na
- K
- Ca
- Zn

---

## JSON/MongoDB input example

Battflow reads electrolyte compositions from documents inside the collections specified in `default.yaml`. A list of JSON files can also be provided.

Battflow primarily reads the electrolyte composition from the SMILES fields. The names provided in `components` are mainly used for visualization of the final results. The exception is the `ions` field, where the user should provide the correct elemental symbol.

```json
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
    "cations": ["0.2M"],
    "ions": ["1M"]
  },
  "simulation_data": {
    "diffusion_coefficients": null,
    "transference_number": null,
    "solvation_statistics": null,
    "coordination_number": null,
    "pairing_percentage": null,
    "dft_energies": null
  }
}
```

---

## Output fields populated by Battflow

### 1. Molecular dynamics transport and solvation properties

The field `simulation_data` will be populated with:

- Diffusion coefficients,
- Transference number,
- The three most common solvation structures,
- Their Cartesian coordinates (XYZ format),
- Coordination number for the solvated ion,
- Pairing percentage for the solvated ion.

---

### 2. DFT energies of individual components

The field `dft_energies` will store:

- DFT energy,
- HOMO,
- LUMO,
- HOMO–LUMO gap

for each molecular component in the electrolyte formulation.
