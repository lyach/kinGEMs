# CASE 1: growth constrained (apply_medium_irreversible(growth_value=avg_growth)), single enzyme (run_optimization_with_dataframe(multi_enzyme_off=False, isoenzymes_off=True,promiscuous_off=True, complexes_off=True)) 
--- Data Preparation (steps 1-3, run once) ---
  Biomass reaction: BIOMASS_Ec_iML1515_core_75p37M
  Baseline growth (unconstrained): 0.8770
=== Step 1: Preparing model data ===
  Found existing files, loading cached data
Warning: Model has 8 validation errors
  Model: 1516 genes, 3682 reactions
=== Step 2: Merging substrate and sequence data ===
  Merged data: 12229 rows
=== Step 3: Processing CPI-Pred kcat values ===
  Processed data: 11789 rows
  Data prep finished in 35.5s
Loaded 3 conditions from davidi_media.csv
## 1/3 Condition: glc_chem_mu_011_v  |  growth=0.11
  Closed 331 reverse-exchange reactions
  Opened 20 exchange reactions from medium specification
  Fixed growth BIOMASS_Ec_iML1515_core_75p37M: lb=ub=0.11
  Baseline growth (media-constrained): 0.1100
  Result: FEASIBLE  |  COBRApy: 0.11  |  EC biomass: 0.11  |  4.3s
## 2/3 Condition: gal_batch_mu_026_s  |  growth=0.26
  Closed 331 reverse-exchange reactions
  Opened 21 exchange reactions from medium specification
  Fixed growth BIOMASS_Ec_iML1515_core_75p37M: lb=ub=0.26
  Baseline growth (media-constrained): 0.2600
  Result: FEASIBLE  |  COBRApy: 0.26  |  EC biomass: 0.26  |  3.4s
## 3/3 Condition: ace_batch_mu_03_s  |  growth=0.3
  Closed 331 reverse-exchange reactions
  Opened 21 exchange reactions from medium specification
  Fixed growth BIOMASS_Ec_iML1515_core_75p37M: lb=ub=0.3
  Baseline growth (media-constrained): 0.3000
  Result: FEASIBLE  |  COBRApy: 0.3  |  EC biomass: 0.3  |  3.9s
=== Pipeline Complete ===
Conditions tested:  3
Feasible:           3/3
# CASE 2: growth unconstrained (apply_medium_irreversible(growth_value=None)), single enzyme (run_optimization_with_dataframe(multi_enzyme_off=False, isoenzymes_off=True,promiscuous_off=True, complexes_off=True)) 
--- Data Preparation (steps 1-3, run once) ---
  Biomass reaction: BIOMASS_Ec_iML1515_core_75p37M
  Baseline growth (unconstrained): 0.8770
=== Step 1: Preparing model data ===
  Found existing files, loading cached data
Warning: Model has 8 validation errors
  Model: 1516 genes, 3682 reactions
=== Step 2: Merging substrate and sequence data ===
  Merged data: 12229 rows
=== Step 3: Processing CPI-Pred kcat values ===
  Processed data: 11789 rows
  Data prep finished in 36.0s
Loaded 3 conditions from davidi_media.csv
## 1/3 Condition: glc_chem_mu_011_v  |  growth=0.11
  Closed 331 reverse-exchange reactions
  Opened 20 exchange reactions from medium specification
  Baseline growth (media-constrained): 26.1161
  Result: FEASIBLE  |  COBRApy: 26.116129  |  EC biomass: 2.721675  |  4.0s
## 2/3 Condition: gal_batch_mu_026_s  |  growth=0.26
  Closed 331 reverse-exchange reactions
  Opened 21 exchange reactions from medium specification
  Baseline growth (media-constrained): 25.7423
  Result: FEASIBLE  |  COBRApy: 25.742258  |  EC biomass: 1.519706  |  3.6s
## 3/3 Condition: ace_batch_mu_03_s  |  growth=0.3
  Closed 331 reverse-exchange reactions
  Opened 21 exchange reactions from medium specification
  Baseline growth (media-constrained): 8.1354
  Result: FEASIBLE  |  COBRApy: 8.135418  |  EC biomass: 1.87319  |  4.7s
=== Pipeline Complete ===
Conditions tested:  3
Feasible:           3/3
# CASE 3: growth unconstrained (apply_medium_irreversible(growth_value=None)), all enzymes (run_optimization_with_dataframe(multi_enzyme_off=False, isoenzymes_off=False,promiscuous_off=False, complexes_off=False))
--- Data Preparation (steps 1-3, run once) ---
  Biomass reaction: BIOMASS_Ec_iML1515_core_75p37M
  Baseline growth (unconstrained): 0.8770
=== Step 1: Preparing model data ===
  Found existing files, loading cached data
Warning: Model has 8 validation errors
  Model: 1516 genes, 3682 reactions
=== Step 2: Merging substrate and sequence data ===
  Merged data: 12229 rows
=== Step 3: Processing CPI-Pred kcat values ===
  Processed data: 11789 rows
  Data prep finished in 34.4s
Loaded 3 conditions from davidi_media.csv
## 1/3 Condition: glc_chem_mu_011_v  |  growth=0.11
  Closed 331 reverse-exchange reactions
  Opened 20 exchange reactions from medium specification
  Baseline growth (media-constrained): 26.1161
  Result: FEASIBLE  |  COBRApy: 26.116129  |  EC biomass: 0.022992  |  4.6s
## 2/3 Condition: gal_batch_mu_026_s  |  growth=0.26
  Closed 331 reverse-exchange reactions
  Opened 21 exchange reactions from medium specification
  Baseline growth (media-constrained): 25.7423
  Result: FEASIBLE  |  COBRApy: 25.742258  |  EC biomass: 0.019438  |  4.6s
## 3/3 Condition: ace_batch_mu_03_s  |  growth=0.3
  Closed 331 reverse-exchange reactions
  Opened 21 exchange reactions from medium specification
  Baseline growth (media-constrained): 8.1354
  Result: INFEASIBLE  |  COBRApy: 8.135418  |  EC biomass: N/A  |  4.4s
=== Pipeline Complete ===
Conditions tested:  3
Feasible:           2/3