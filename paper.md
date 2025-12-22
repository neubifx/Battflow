---
title: 'Battflow: an automated workflow for properties prediction of liquid electrolytes'
tags:
  - Python
  - batteries
  - electrolyte
  - DFT
  - molecular dynamics
authors:
  - name: Neubi F. Xavier Jr
    orcid: 0000-0002-2133-0557
    corresponding: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Matthias J. Gollomb
    orcid: 0000-0001-6749-0129
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Qiong Cai
    orcid: 0000-0002-1677-0515
    corresponding: true
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: School of Chemistry and Chemical Engineering, University of Surrey, Guildford, GU2 7XH, United Kingdom
   index: 1
   #ror: 00hx57361

date: 22 December 2025
bibliography: paper.bib
---

# Summary

Battflow provide an automated workflow, integrating a suite of Python packages, GROMACS [@abraham_2025_17671776]
and ORCA [@ORCA] to predict properties of liquid electrolytes for Li metal batteries. It aims to simplify the user 
interactions when generating large MD/DFT datasets, aiding the high-throughput screening and discovery of new electrolytes.
The workflow is linked to collection inside a MongoDB cluster in which missing calculated properties are identified and 
flagged for further calculations. The only information needed to start the simulations are the SMILES string and molar concentration 
of each electrolyte component. A document for each electrolyte composition is then created and diffusion coefficients and solvation 
structure statistics are computed from molecular dynamics simulations and binding energies for solvation clusters and HOMO-LUMO 
energies of each individual components are updated. 

# Statement of need

Lithium metal batteries (LMB) are regarded as a promising solution to meet the market demand for energy storage systems with 
high specific capacity. However, the implementation of lithium metal anodes is currently hindered by poor cycle life and 
uncontrollable side reactions between Li metal and liquid electrolytes. Liquid electrolyte engineering, which involves mixing 
different molecules to create electrolytes with specific properties, is ultimately the most cost-effective approach for making
LMBs viable. However, there are a myriad of possible electrolyte formulations due to the large number of commercially available
molecules, recently synthesized electrolyte-specific compounds, and various strategies for fine-tuning electrolyte components.
Consequently, there is an increased need for theory-guided rational design of new electrolyte formulations for LMBs, aiming to 
reduce research costs and avoid “trial-and-error” approaches.

There is an increased need of standardized computational data that can guide the experiments of liquid electrolytes for lithium 
metal batteries. A combination of molecular dynamics and Density Functional Theory (DFT) simulations for the estimation of transport
and electronic properties of the bulk electrolyte and their individual components shows the best balance between accuracy and cost
efficiency for proeprties prediction. However, setting up force-fields and settings for molecular dynamics simulations, followed by 
DFT calculations of relevant Li solvation clusters, for several electrolyte composition can be a daunting task. Battflow is intended to 
both theoreticians and experimentalists to be able to install it on a Linux environment, already providing a out-of-the-box default 
settings to run the workflow with basic setttings.

# Usage and availability

Battflow inputs are documents inside a collection on MongoDB. The information to connect to MongoDB, either through localhost or 
is provided inside `config.yaml`. A `.json` example input which should be uploaded in MongoDB as a document is provided as an example.
Battflow reads from `smiles` and `concentrations` fields inside each document to build the molecular structure for each component and to 
create the electrolyte box, respectively. The output is reported in the `simulation_data` field, in which if any of the calculated proeprty 
is flagged as absent, the workflow is flagged to start. Molecular dynamics steps are 1) The setup of the GAFF2 force-fields for the molecular, 
created by ACPYPE [@SousadaSilva2012], 2) Creation of the electrolyte box 3) running of the simulations following a minimization, equilibration 
and production runs using GROMACS [@abraham_2025_17671776] and 4) uploading of the diffusion properties and solvations stucture statistics into 
the MongoDB document. The DFT simulations in ORCA [@ORCA] follows the output of the molecular dynamics runs, in which calculations are done
respective to the 3 more present solvation clusters. Binding energy of each cluster, as well as HOMO-LUMO energies for each component is computed and
uploaded back to MongoDB. 

Battflow is available for Linux operating systems and can be downloaded from GitHub (https://github.com/neubifx/Battflow/tree/main) through the 
GPL-3.0 license. Additional documentations is already available in the repository page and is being constantly updated.

# Acknowledgements

The work was funded by the AIchemy Pump-Priming project (EPSRC\\XXXXXX). The authors acknowledge funding by Horizon Europe through the OPERA consortium
(Grants Number 101103834) and under the UKRI Horizon Europe Guarantee Extension (Ref Number 10078555), by the Faraday Institution through the LiSTAR program 
(Grants FIRG014, FIRG058), and by Royal Society (IEC\\NSFC\\211200).

# References
