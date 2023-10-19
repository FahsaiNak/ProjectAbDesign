# ProjectAbDesign

## About
This project transforms the antibody design framework developed by Aguilar Rangel et al. [cite](https://doi.org/10.1126/sciadv.abp9540) into an automated workflow.

## Getting started

### Prerequisites

* Install [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) with python 3.11

* Install [Master](https://grigoryanlab.org/master/): a rapid structural similarity search program

### Installation

1. Clone the repository
   ```sh
   git clone https://github.com/FahsaiNak/ProjectAbDesign.git
   ```

2. Create environment from environment.yml
   ```sh
   conda env create -f environment.yml
   ```

### Step 1: Antibody-Antigen database construction
Antibody-Antigen (AbAg) database is a collection of antigen-like and CDR-like regions in all non-redundant general protein structures reported in the [PDB database](https://www.rcsb.org/docs/programmatic-access/file-download-services).

1. PDB90 database

2. CDR database and fragmented-CDR datasets

3. CDR-like region identification

4. Antigen-like region identification

5. AbAg Database

### Step 2: Antigen epitope processing and its interacting-CDR-like searching

## Workflow

### 1. AbAg Database construction (Fahsai, David, Juan)
  -	PDB90 dataset (PDB), CDR-fragment dataset (full Ab, SabDab)
  -	Blastp + MASTER for searching CDR-like regions in the PDB90 dataset (target) by CDR-fragment dataset (query)
  -	CACA + SASA to find the interacting (Ag-like) regions of the CDR-like regions within the same protein target
  -	Store the database

### 2. CDR-like identification (David, Lindsey)
  - **Input = full epitope structure, Output = CDR-like candidate structures**
  -	Epitope fragments (target antigen), AbAg database
  -	Blastp + MASTER for searching for Ag-like regions in the AbAg database (target) by Epitope fragments (query)
  -	Track back to the corresponding CDR-like regions of the matching Ag-like regions

### 3. CDR selection and optimization (Lindsey, Juan)
  -	**Input = CDR-like candidate structures, Output = Optimized CDR candidate structures**
  -	Combine, rank and select CDR-like regions by their interactions with the epitopes and some properties (solubility)
  -	(optional) Contact residual optimization

### 4. CDR grafting (Fahsai, ________)
  -	**Input = Optimized CDR candidate structures, Output = CDR-grafted nanobodies (hits)**
  -	Optimized CDRs, nanobody-CDR dataset (SabDab)
  -	MASTER for searching for structurally-matching nanobody-CDRs in nanobody-CDR dataset (CDR1 and/or CDR3, target) by the optimized CDRs (query)
  -	Replace the matching regions with the optimized CDRs with all possible combinations
  -	Contact residual optimization
  -	Rank and select the hits
