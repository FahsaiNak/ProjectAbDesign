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
- Ensure you have enough free space for the PDB90 files of interest.
- Add a CSV list of target PDB proteins to download in the Datasets directory. The list most have a similar format to the somePDB.csv file in the Datasets directory.
- Once the sequence list has been added, modify the PDB90.sh script in the run folder. It currently has python ../src/PDB90.py --output_folder "../Datasets/all_PDB" --csv_file "../Datasets/somePDB.csv" written, change "../Datasets/somePDB.csv" with the location of the downloaded PDB csv file.
- Move to the run directory and run the bash script with:
  ```sh
  bash PDB90.sh
  ```
- The script will take a long time to run, as there are around 40000 .pdb files that are being downloaded, uncompressed and cleaned.

2. CDR database and fragmented-CDR datasets
- Navigate to [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/search/?all=true#downloads)
- Click on Downloads on the left hand side of the page
- Click on Download an archived zip file to download structures
- Once the zip file is downloaded, extract the chothia zip file and save it into the Datasets directory in the local repository
- Change directory to the run directory
- On the command line run:
  ```sh
  bash CDR_fragment_database.sh
  ```

3. CDR-like region identification
  - Remain in run directory
  - Run Prep_queries_AbAg.sh and Prep_targets_AbAg.sh to convert fragmented-CDRs (queries) and cleaned-PDB90 proteins (targets) in PDB to PDS that is readable for Master program.
    ```sh
    ./Prep_queries_AbAg.sh
    ./Prep_targets_AbAg.sh
    ```
  - Run Ablike.sh search for antigen(CDR)-like regions in every target proteins (PDB90 proteins)
    ```sh
    ./Ablike.sh
    ```
  - Run get_matchInfo.sh to collect antigen(CDR)-like structure information from the Master match files.
    ```sh
    ./get_matchInfo.sh
    ```

4. Antigen-like region identification and AbAg database generation
  - Remain in run directory
  - Run Aglike.sh to search for contacting/interacting residues of each antigen(CDR)-like region in the same protein. Then, all the antibody-antigen-like information from each PDB90 proteins is combined into a compressed file named AbAb.pkl
    ```sh
    ./Aglike.sh
    ```

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
