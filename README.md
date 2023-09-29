# ProjectAbDesign

## Workflow

### AbAg Database construction (Fahsai, ________, ________)
  -	PDB90 dataset (PDB), CDR-fragment dataset (full Ab, SabDab)
  -	Blastp + MASTER for searching CDR-like regions in the PDB90 dataset (target) by CDR-fragment dataset (query)
  -	CACA + SASA to find the interacting (Ag-like) regions of the CDR-like regions within the same protein target
  -	Store the database

### CDR-like identification (________, ________)
  - **Input = full epitope structure, Output = CDR-like candidate structures**
  -	Epitope fragments (target antigen), AbAg database
  -	Blastp + MASTER for searching for Ag-like regions in the AbAg database (target) by Epitope fragments (query)
  -	Track back to the corresponding CDR-like regions of the matching Ag-like regions

### CDR selection and optimization (________, ________)
  -	**Input = CDR-like candidate structures, Output = Optimized CDR candidate structures**
  -	Combine, rank and select CDR-like regions by their interactions with the epitopes and some properties (solubility)
  -	(optional) Contact residual optimization

### CDR grafting (Fahsai, ________)
  -	**Input = Optimized CDR candidate structures, Output = CDR-grafted nanobodies (hits)**
  -	Optimized CDRs, nanobody-CDR dataset (SabDab)
  -	MASTER for searching for structurally-matching nanobody-CDRs in nanobody-CDR dataset (CDR1 and/or CDR3, target) by the optimized CDRs (query)
  -	Replace the matching regions with the optimized CDRs with all possible combinations
  -	Contact residual optimization
  -	Rank and select the hits
