
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![License](https://img.shields.io/badge/license-GPLv3-green)
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

Battflow relies on [ACPYPE](https://github.com/alanwilter/acpype) for GAFF2 topology generation. A Linux-based system is recommended for installation and execution.

Clone the repository:

```bash
git clone https://github.com/neubifx/Battflow.git
cd Battflow
```

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate battflow_env
```

Install Battflow:

```bash
pip install .
```

For development purposes, Battflow can alternatively be installed in editable mode:

```bash
pip install -e .
```

---

### 2. External software requirements

Battflow requires the following external software:

- [MongoDB](https://www.mongodb.com/)
- [GROMACS](https://www.gromacs.org/)
- [ORCA](https://orcaforum.kofo.mpg.de/)
- [AmberTools](https://ambermd.org/AmberTools.php) (`antechamber` and `tleap`)
- [Packmol](https://m3g.github.io/packmol/)

These programs must be installed and accessible in the execution environment. The required Python packages, including ACPYPE, are installed through the provided environment and package configuration.

Ensure that the GROMACS command, AmberTools executables and Packmol are available in the system environment, and that the ORCA executable path is correctly configured.

---

## Configuration

Users can inspect and modify the default configuration according to their simulation environment. Generate a copy of the default configuration file in the current working directory:

```bash
battflow --write-config
```

This creates:

```text
custom_battflow.yaml
```

Users can then modify this file according to their local environment (e.g., MongoDB URI, ORCA and GROMACS commands, HPC submission commands, number of cores, etc.).

---

### MongoDB configuration

This section should point to the database and collection from which Battflow will retrieve electrolyte compositions for calculation. The default options correspond to a local MongoDB installation, although a MongoDB Atlas URI can also be provided in `host`.

```yaml
mongodb:
    host: localhost
    port: 27017
    database: working_db
    collection: working_collection
```

---

### MD run environment

The default options assume a standard Linux installation of GROMACS. These commands can be adapted to match the commands used in an HPC environment.

```yaml
md_run_env:
    mdrun: "gmx mdrun"   # Command used to run GROMACS
    mpiexec: ""          # MPI executor, if required (e.g., srun)
    ncores: 8            # Number of MPI processes when using an MPI executor
```

When no MPI executor is specified, GROMACS determines the number of threads according to its own settings.

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
battflow --config custom_battflow.yaml
```

Optional: specify the solute ion used during MD analysis. If no option is provided, Li will be considered the solute ion.

```bash
battflow --config custom_battflow.yaml --solute-ion li
```

Supported solute ions currently include:

- Li
- Na
- K
- Ca
- Zn

---

## JSON/MongoDB input example

Battflow reads electrolyte compositions from documents inside the collection specified in `custom_battflow.yaml`. JSON documents can be imported into the MongoDB collection before running the workflow.

Battflow primarily reads the electrolyte composition from the SMILES fields. The names provided in `components` are mainly used for visualization of the final results. The exception is the `ions` field, where the user should provide the correct elemental symbol.

The following example contains an electrolyte composition with ethylene carbonate (EC), ethyl methyl carbonate (EMC), hexafluorophosphate (PF6), methylimidazolium and lithium.

```json
{
  "components": {
    "molecules": ["ec", "emc"],
    "anions": ["pf6"],
    "cations": ["methylimidazolium"],
    "ions": ["li"]
  },
  "smiles": {
    "molecules": ["C1COC(=O)O1", "COC(=O)OCC"],
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

The `simulation_data` fields are initially set to `null` to indicate properties that have not yet been calculated.

MongoDB can automatically generate the document's `_id` during import.

---

## Output fields populated by Battflow

### 1. Molecular dynamics transport and solvation properties

The field `simulation_data` will be populated with:

- Diffusion coefficients
- Transference number
- The three most common solvation structures
- Their Cartesian coordinates (XYZ format)
- Coordination number for the solvated ion
- Pairing percentage for the solvated ion

---

### 2. DFT energies of individual components

The field `simulation_data.dft_energies` will store:

- DFT energy
- HOMO
- LUMO
- HOMO–LUMO gap

for each molecular component in the electrolyte formulation.
