---
title: 'Battflow: an automated workflow for predicting key properties of battery electrolytes'
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
  - name: Matthias J. Golomb
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

date: 11 May 2026
bibliography: paper.bib
---

# Summary

`Battflow` provides an automated workflow that integrates a suite of Python packages, GROMACS [@abraham_2025_17671776], and ORCA [@ORCA] to predict properties of battery electrolytes. The workflow aims to simplify the generation of large MD/DFT datasets, thereby supporting the high-throughput screening and discovery of new electrolyte formulations. `Battflow` is directly connected to a MongoDB collection, where missing calculated properties are automatically identified and flagged for further calculations. The only information required to initiate the simulations is the SMILES string and molar concentration of each electrolyte component. After identifying electrolyte compositions with missing calculated properties from MongoDB documents, `Battflow` automatically generates the electrolyte simulation box and performs molecular dynamics simulations to compute diffusion coefficients and solvation structure statistics. Subsequently, density functional theory calculations are performed to estimate HOMO–LUMO energies of the individual electrolyte components. The calculated properties are then automatically uploaded back to MongoDB.

# Statement of need

Lithium metal batteries (LMBs) are regarded as a promising solution for meeting the growing demand for energy storage systems with high specific capacity. Several battery technologies incorporating lithium metal anodes have recently attracted significant attention, including lithium-sulfur (Li-S), lithium-oxygen (Li-O~2~), and lithium-carbon dioxide (Li-CO~2~) batteries. However, the implementation of lithium metal anodes is currently hindered by poor cycle life and uncontrolled side reactions between Li metal and liquid electrolytes. Liquid electrolyte engineering, which involves combining different molecular components to create electrolytes with tailored properties, is considered one of the most cost-effective strategies for enabling commercially viable LMBs. Nevertheless, the number of possible electrolyte formulations is extremely large due to the wide range of commercially available molecules, newly synthesised electrolyte-specific compounds, and multiple strategies for fine-tuning electrolyte compositions. Consequently, there is an increasing need for theory-guided rational design of new electrolyte formulations for LMBs to reduce research costs and minimise trial-and-error approaches.

There is also an increasing need for standardised computational data to support experimental studies on battery electrolytes, particularly in light of the rapid growth of the literature in this field, with approximately 200 papers containing the keywords “battery electrolyte” published per week in 2023 alone. A combination of molecular dynamics (MD) and density functional theory (DFT) simulations for the estimation of transport and electronic properties of both bulk electrolytes and their individual components offers a good balance between predictive accuracy and computational cost. However, setting up force fields and simulation parameters for molecular dynamics simulations, followed by DFT calculations of relevant Li solvation clusters across multiple electrolyte compositions, can be a daunting and error-prone task even for experienced researchers.

`Battflow` is intended for both theoreticians and experimentalists with varying levels of expertise in atomistic modelling tools and provides out-of-the-box default settings to execute automated workflows using only basic configuration inputs. At present, `Battflow` extends beyond Li-metal and Li-ion batteries and also supports simulations of Na-, K-, Zn- and Ca-based battery electrolytes.

# State of the field

Obtaining curated data is particularly challenging in rapidly evolving research fields due to the continuous growth in the number of published studies. In 2024 alone, approximately 240,000 publications included keywords related to batteries, making it impossible for a single research group to comprehensively follow all newly released information [@CHEN2026101091]. Consequently, there is an increasing need for curated and easily accessible data, following FAIR (Findability, Accessibility, Interoperability, and Reuse) principles, to improve the understanding of battery performance across a wide range of operating conditions. Although significant advances have been made in high-throughput experimental techniques aimed at accelerating the generation of battery materials data, important limitations remain regarding the diversity of materials that can be synthesised and the range of properties that can be experimentally measured [@Xu31122024]. Computationally derived data can therefore play an important role in supplementing missing information while reducing the trial-and-error nature of experimental exploration.

However, computational workflows used to predict battery electrolyte properties are typically operated by specialists in computational modelling. In many cases, molecular dynamics (MD) simulations and density functional theory (DFT) calculations are performed independently using different software packages. MD simulations are commonly conducted using open-source software such as GROMACS or LAMMPS, while DFT calculations are performed using packages such as Gaussian or ORCA. These approaches rely on different methodologies and require distinct expertise, making it challenging for users familiar with one methodology to seamlessly adopt the other. This barrier is even greater for experimental researchers seeking to predict battery electrolyte properties. As a result, computational data reported in the literature is often limited to isolated properties, such as HOMO–LUMO energies or MD-derived transport properties, providing only partial insight into electrolyte behaviour.

`Battflow` addresses this gap by providing an accessible workflow that integrates MD and DFT simulations within a single framework for battery electrolyte systems. The workflow enables out-of-the-box usage requiring only basic information about the electrolyte composition, lowering the barrier for experimental researchers and newcomers to computational modelling. At the same time, advanced users retain full control over simulation parameters through configurable inputs, allowing the workflow to be adapted for specialised or high-throughput studies.

# Software design

`Battflow` is designed to provide out-of-the-box usability while remaining flexible for deployment on high-performance computing (HPC) systems. The workflow targets large-scale electrolyte simulations and therefore prioritises integration with established open-source simulation engines: GROMACS for classical molecular dynamics simulations and ORCA for density functional theory calculations, reflecting the common requirement for large numbers of processors when modelling complex electrolyte formulations.

GROMACS simulations are executed through a dedicated wrapper (`GromacsWrapper`), while ORCA calculations are orchestrated through the `Atomic Simulation Environment (ASE)` calculator interface. A modified class adapts job submission to different HPC schedulers and execution environments according to user-defined settings specified through the `battflow --write-config` command.

Within the generated configuration file, `Battflow` default settings can be modified to adapt the workflow to different computational environments and software installations, including different MongoDB `host`, `database`, and `collection` settings (including MongoDB Atlas), execution commands for GROMACS and ORCA, levels of theory adopted in DFT calculations, and the number of processors used during MD and DFT simulations. `Battflow` can then be executed with user-defined settings using `battflow --config custom_battflow.yaml`.

To enable computationally ready data management, `Battflow` is designed for direct integration with MongoDB, supporting both local and cloud-based installations depending on internet accessibility on computational nodes. The workflow is orchestrated entirely in Python within an Anaconda environment, allowing straightforward installation of supporting tools such as ACPYPE for AMBER force-field generation. Python-based analysis libraries, including MDAnalysis and SolvationAnalysis, are used for post-processing and data extraction.

# Research impact statement

`Battflow` has already been actively used despite its relatively recent development. It is currently the primary workflow employed for the computational data generation of nearly 100 electrolyte formulations as part of a combined experimental–computational dataset extracted from the literature using large language models (LLMs), which is close to public release. In addition, one scientific publication in which `Battflow` was used to simulate Ca-ion battery electrolytes has already been published [@Luo2026], and a preprint is available online [@Soni2025].

The package has additionally attracted contributions from a developer beyond the main author (\@neubifx). Beyond research applications, components of the workflow have been adapted into Jupyter notebooks for undergraduate Chemistry teaching, enabling students to construct molecular dynamics simulation boxes, execute simulations, perform analysis, and visualise results within a single environment.

Together, these applications demonstrate that `Battflow` serves both as an accessible entry point for experimentalists and early-career researchers new to computational methods, and as a scalable tool for high-throughput data generation in battery electrolyte research.

# Usage and availability

`Battflow` is pip-installable and is configured by default to support standard Linux installations of GROMACS, ORCA, and local MongoDB instances. For quick dependency and installation using conda, run the following commands:

`conda env create -f environment.yml`
`conda activate battflow_env`

`pip install git+https://github.com/neubifx/Battflow.git`

A modified class adapts job submission to different HPC schedulers and execution environments according to the settings specified in the `custom_battflow.yaml` configuration file, which can be automatically generated and edited using `battflow --write-config`. This command creates a user-editable configuration file (`custom_battflow.yaml`) in the working directory. `Battflow` can then be executed using:

`battflow --config custom_battflow.yaml`

A `.json` example input file, which should be uploaded to MongoDB as a document, is provided in the repository. `Battflow` reads the `smiles` and `concentrations` fields within each document to build the molecular structures of the electrolyte components and generate the electrolyte box, respectively. The calculated properties are stored in the `simulation_data` field. If any target property in a MongoDB document is flagged as absent, the workflow automatically triggers the corresponding simulations and analyses.

The molecular dynamics workflow consists of: (1) generation of GAFF2 force fields for the molecular components through ACPYPE [@SousadaSilva2012]; (2) construction of the electrolyte box; (3) execution of minimisation, equilibration, and production molecular dynamics runs using GROMACS [@abraham_2025_17671776]; and (4) extraction and uploading of diffusion properties and solvation structure statistics to MongoDB. DFT simulations using ORCA [@ORCA] are subsequently performed for the three most populated solvation clusters identified during the molecular dynamics simulations. HOMO, LUMO, and HOMO–LUMO gap energies are calculated for both the complete solvation clusters and their isolated molecular components. After the analysis is completed, the results are uploaded back to MongoDB.

`Battflow` is available for under the GPL-3.0 licence. Additional documentation is available in the repository and is continuously updated.

# AI usage disclosure

Generative AI tools were used to assist with the initial drafting of descriptions for selected functions and to perform grammar and style corrections in the documentation. All AI-generated text and function descriptions were thoroughly reviewed and verified by the authors prior to submission.

# Acknowledgements

This project has received funding from the AI for Chemistry: AIchemy Hub (EPSRC grant EP/Y028775/1 and EP/Y028759/1). The authors acknowledge funding from Horizon Europe through the OPERA consortium (Grant Number 101103834) and under the UKRI Horizon Europe Guarantee Extension (Ref Number 10078555), from the Faraday Institution through the LiSTAR programme (Grants FIRG014, FIRG058), and from the Royal Society (IEC\\NSFC\\211200).

# References
