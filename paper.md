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

date: 22 December 2025
bibliography: paper.bib
---

# Summary

`Battflow` provides an automated workflow that integrates a suite of Python packages, GROMACS [@abraham_2025_17671776], and ORCA [@ORCA] to predict properties of battery electrolytes. It aims to simplify user interactions when generating large MD/DFT datasets, thereby aiding the high-throughput screening and discovery of new electrolytes. The workflow is linked to a collection within a MongoDB database, in which missing calculated properties are identified and flagged for further calculations. The only information needed to start the simulations is the SMILES string and molar concentration of each electrolyte component. After identifying the electrolyte composition with missing calculated properties from MongoDB documents, an electrolyte box is then created, and diffusion coefficients and solvation structure statistics are computed from molecular dynamics simulations, while binding energies for solvation clusters and HOMO–LUMO energies of each individual component are updated.

# Statement of need

Lithium metal batteries (LMBs) are regarded as a promising solution for meeting market demand for energy storage systems with high specific capacity. Recently, several battery technologies incorporating lithium metal anodes have attracted increasing attention, including lithium-sulfur (Li-S) batteries, lithium-oxygen (Li-O~2~), and lithium-carbon dioxide (Li-CO~2~) batteries. However, the implementation of lithium metal anodes is currently hindered by poor cycle life and uncontrollable side reactions between Li metal and liquid electrolytes. Liquid electrolyte engineering, which involves mixing different molecules to create electrolytes with specific properties, is ultimately the most cost-effective approach for making LMBs viable. However, there are myriad possible electrolyte formulations due to the large number of commercially available molecules, recently synthesised electrolyte-specific compounds, and various strategies for fine-tuning electrolyte components. Consequently, there is an increased need for theory-guided rational design of new electrolyte formulations for LMBs, aiming to reduce research costs and avoid “trial-and-error” approaches.

There is an increasing need for standardised computational data to guide experimental studies on battery electrolytes, particularly in light of the rapid growth of the literature in this field, with approximately 200 papers containing the keywords “battery electrolyte” published per week in 2023 alone. A combination of molecular dynamics (MD) and Density Functional Theory (DFT) simulations for the estimation of transport and electronic properties of the bulk electrolyte and its individual components offers a good balance between accuracy and cost efficiency for properties prediction. However, setting up force fields and simulation settings for molecular dynamics, followed by DFT calculations of relevant Li solvation clusters for several electrolyte compositions, can be a daunting task and error-prone task, even for experienced researchers. 

`Battflow` is intended for both theoreticians and experimentalists with different levels of expertise using atomistic modelling tools and can provide out-of-the-box default settings to run the automated workflow with basic configurations. At present, `Battflow` extends beyond Li-metal and Li-ion batteries and supports the usage of other battery systems, including Na, K, Zn, Ca and Mg-based batteries.

# State of the field  

Obtaining curated data is particularly challenging in “hot” research fields due to the rapid growth in the number of published papers. In 2024 alone, approximately 240,000 publications included keywords related to batteries, making it impossible for a single research group to follow all the information being released [@XCHEN2026101091]. Therefore, there is an increase need on curated and easilly accesible data, accordingly to FAIR (Findability, Accessibility, Interoperability, and Reuse) standards to fully understand battery performance at a wide range of operation conditions. Although significant advancements have been made in high-throughput experimental techniques aimed at accelerating the generation of battery materials data, challenges remain regarding the diversity of materials that can be synthesised and the range of properties that can be measured [@Xu31122024]. Computationally derived data can therefore play an important role in supplementing missing information while reducing the trial-and-error nature of experimental exploration.

However, computational workflows used to predict battery electrolyte properties are typically operated by specialists in computational modelling. In many cases, MD simulations and DFT calculations are performed independently using different software packages. MD simulations are commonly conducted using open-source software such as GROMACS or LAMMPS, while DFT calculations are performed using packages such as GAUSSIAN or ORCA. These approaches rely on different methodologies and require distinct expertise, making it challenging for users of one methodology to adopt the other seamlessly. This barrier is even higher for experimental researchers seeking to predict battery electrolyte properties. As a result, computational data reported in the literature is often limited to isolated properties, such as HOMO–LUMO energies or MD-derived transport properties, providing only partial insight into electrolyte behaviour.

``Battflow`` addresses this gap by providing an accessible workflow that integrates MD and DFT simulations within a single framework for battery electrolyte systems. The workflow enables out-of-the-box usage requiring only basic information about the electrolyte composition, lowering the barrier for experimental researchers and newcomers to computational modelling. At the same time, advanced users retain full control over simulation parameters through configurable inputs, allowing the workflow to be adapted for more specialised or high-throughput studies.

# Software design

`Battflow` is designed to provide out-of-the-box usability while remaining flexible for deployment on high-performance computing (HPC) systems. The workflow targets large-scale electrolyte simulations and therefore prioritises integration with established open-source simulation engines: GROMACS for classical molecular dynamics and ORCA for density functional theory calculations, reflecting the common requirement for large numbers of processors when modelling complex electrolyte formulations. `Battflow` installation can be pip installable and it will be setted to fit default installation settings of GROMACS and ORCA and a local installation of MongoDB:

`conda env create -f environment.yml
conda activate battflow_env`

`pip install git+https://github.com/neubifx/Battflow.git`


GROMACS simulations are executed through a dedicated wrapper (GromacsWrapper) and ORCA calculations are orchestrated via the Atomic Simulation Environment (ASE) calculator interface. A modified class adapts job submission to different HPC schedulers and execution environments as specified by configuration in the `default.yaml` file which can be simple modified by calling `battflow --write-config`. Within this file which is copied in the main folder, `Battflow` default settings can be changed to adapt to different environments and installations, e.g. different `host`, `database` and `collections` in MongoDB (including Atlas), execution command of GROMACS and ORCA, level of theory adopted in DFT simulations and number of processors adopted in the DFT and MD simulatons. Then, `Battflow` can be simply run with custom configuration settings by using:

`battflow --config custom_battflow.yaml`

To enable computationally ready data management, `Battflow` is designed for direct integration with MongoDB, supporting both local and cloud-based installations depending on the internet accessibility on the computational nodes. The workflow is orchestrated entirely in Python within an Anaconda environment, allowing straightforward installation of supporting tools such as ACPYPE for AMBER force-field generation. Python-based analysis libraries, including MDAnalysis and Solvation Structure, are used for post-processing and data extraction. Results are uploaded back in the document containing the electrolyte information in MongoDB, including `xyz` coordinates with the most present solvation structures in the elecrolyte, guaranteeing standardization and improved readability of the data generated.

# Research impact statement

`Battflow` has been actively used despite its relatively recent development. It is currently the primary workflow employed for the computational data generation of nearly 100 electrolyte formulations as part of a combined experimental–computational dataset extracted from the literature using LLMs, which is close to public release. In addition, one scientific publications in which `Battflow` was used to perform simulations of alkali-metal battery electrolytes have been published [@Soni2025]. The package has attracted contributions from an additional developer beyond the main author (@neubifx). Beyond research applications, components of the workflow have been adapted into Jupyter notebooks for undergraduate teaching, enabling students to construct molecular dynamics simulation boxes, run simulations, perform analysis, and visualise results within a single environment. Together, these applications demonstrate that `Battflow` serves both as an accessible entry point for experimentalists and early-career researchers new to computational methods, and as a scalable tool for high-throughput data generation in battery electrolyte research.

# Usage and availability

`Battflow` inputs consist of documents stored within a MongoDB collection. The information required to connect to MongoDB, either through localhost or a remote instance, is provided in `config.yaml`. A `.json` example input file, which should be uploaded to MongoDB as a document, is provided. `Battflow` reads the `smiles` and `concentrations` fields within each document to build the molecular structure of each component and to create the electrolyte box, respectively. The output is reported in the `simulation_data` field; if any calculated property is flagged as absent, the workflow is triggered to start.

The molecular dynamics workflow consists of: (1) setup of GAFF2 force fields for the molecular components, generated using ACPYPE [@SousadaSilva2012]; (2) creation of the electrolyte box; (3) execution of minimisation, equilibration, and production MD runs using GROMACS [@abraham_2025_17671776]; and (4) uploading diffusion properties and solvation structure statistics to the MongoDB document. The DFT simulations in ORCA [@ORCA] follow the molecular dynamics runs, with calculations performed for the three most prevalent solvation clusters. Binding energies are computed for every solvation cluster and HOMO-LUMO energies are calculated for the entire cluster and separated components, as well. After the analysis, results are uploaded back to MongoDB.

`Battflow` is available for Linux operating systems and can be downloaded from GitHub (https://github.com/neubifx/Battflow/tree/main) under the GPL-3.0 licence. Additional documentation is available on the repository page and is continuously updated.

# AI usage disclosure

Generative AI tools were used to assist with the initial drafting of descriptions for selected functions and to perform grammar and style corrections in the documentation. All AI-generated text and function descriptions were thoroughly reviewed and verified by the authors prior to submission.

# Acknowledgements

This project has received funding from the AI for Chemistry: AIchemy Hub (EPSRC grant EP/Y028775/1 and EP/Y028759/1). The authors acknowledge funding from Horizon Europe through the OPERA consortium (Grant Number 101103834) and under the UKRI Horizon Europe Guarantee Extension (Ref Number 10078555), from the Faraday Institution through the LiSTAR programme (Grants FIRG014, FIRG058), and from the Royal Society (IEC\\NSFC\\211200).

# References
