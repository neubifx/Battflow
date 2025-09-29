# Battflow - Automated Workflow for Properties Prediction of Liquid Electrolytes

This project aims to automatically identify and calculate missing properties in a MongoDB-based liquid electrolyte database. FireWorks is used to orchestrate the workflow across computational resources.

## Overview

The workflow is composed of the following steps:
1. **Query MongoDB** for entries with missing properties.
2. **Preprocess** structures from SMILES information and setup force-fields and electrolyte composition.
3. **Run MD simulations** using GROMACS.
4. **Perform molecular calculations** using ORCA.
5. **Postprocess** results and update MongoDB.

Each step is implemented as a FireTask and combined into a FireWorks workflow.

## Installation

> _To be completed_

- Python ≥ 3.9
- [MongoDB](https://www.mongodb.com/)
- [ASE](https://wiki.fysik.dtu.dk/ase/)
- [GROMACS](https://www.gromacs.org/)
- [ORCA](https://orcaforum.kofo.mpg.de/)
- [FireWorks](https://materialsproject.github.io/fireworks/)

