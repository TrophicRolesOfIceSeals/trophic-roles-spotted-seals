# trophic-roles-spotted-seals

Welcome to a Github repository associated with the paper, 
"Trophic roles of spotted seals in the Bering, Chukchi, and Beaufort Seas" by
Boveng et al.

Here, we provide data and code to recreate all analyses in this paper.  Data
are provided in the "data" directory, while R markdown and R scripts are located 
in the "inst" directory.

Analysis proceeded in a number of steps, with results from some steps feeding
into others.  It will thus be necessary to run scripts in order to obtain intermediate 
products before proceeding to the final analysis.  It may also be necessary to 
edit paths to various objects in the scripts to ensure they are accessible in 
your local computing environment.

The following is a possible ordering of scripts to recreate the analysis, along with
a description of each.  In some cases, rendered versions of .Rmd files are also provided 
as supplements of the manuscript.  Supplement numbers will differ from the following step
numbers because supplements match the order in which material is discussed in the manuscript,
which may differ from the order in which analysis needs to occur.  For ease of implementation, 
we include *most* of the intermediate products as .RData files in the `data` directory.  However,
owing to the large size of the .RData file associated with seal weight-at-age relationships (which
exceeds github size limits), users must at least run step 6 (supplement 5) before running the 
final integrated analysis.


1. Run `otolith_length_modeling_and_lit.Rmd`, which describes relevant literature and selects
   a preferred otolith length to fish length conversion for each priority fish prey species.
   
    **Inputs**: `/data/otolith_data_for_Paul.csv`
    
    **Outputs**: `/data/otolith_regressions.RData`   

2. Run `mass_energy.Rmd` to produce relationships between fish wet mass and caloric content.

    **Inputs**: `/data/Lt-Wt-EnergyWet.csv`
    
    **Outputs**: `/data/energy_mass_regressions.RData`

3. Run `Len_Wt.R` to produce length-weight relationships for various spotted seal fish prey species.

    **Inputs**: `/data/Lth-Wt-Energy-Combined.csv`
    
    **Outputs**: `/data/weight_length_regressions.RData`

4. [Supplement 3] Run `analyze_stomachs.Rmd` to produce posterior predictions of otolith length and frequency
    found in spotted seal stomach and scat samples.

    **Inputs**: `/data/ADFG_SpottedSeals_OtolithsMeasured20220627.csv`, `/data/AFSC_spottedprey_wotos_meas.csv`,
    `/data/ADFG_SpottedSeals_OtolithsMeasured20230424_PriorityMeasurements.csv`, 
    `/data/ADFG_SpottedSeals_StomachContentsData20220705.csv`
    
    **Outputs**: `/data/otolith_posteriors.RData`, `/data/Ot_len_df.RData`

5. Run `calculate_energy_composition.Rmd`, which uses posterior predictions for the length 
    and expected minimum number of fish represented by otoliths in spotted seal stomachs, 
    otolith length-fish length conversions, fish morphometrics, and energy content conversions to 
    calculate posterior distributions for prey composition.
    
    **Inputs**: `/data/otolith_posteriors.RData`, `/data/Ot_len_df.RData`, `/data/weight_length_regressions.RData`,
    `/data/energy_mass_regressions.RData`, `/data/otolith_regressions.RData`
    
    **Outputs**: `/data/prey_composition_g_kcal.RData`

6. [Supplement 5] Run `seal_mass_at_age.Rmd` to predict weight-at-age of spotted seals using a Bayesian 
    implementation of von Bertalanffy growth models.
    
    **Inputs**: `data/Spotted seal live capture measurements.csv`, `data/Stable_stage_matrix.RData`
    
    **Outputs**: `data/Seal_wt_age.RData`

7.  [Supplement 6] Run `captiveFoodIntake_massSpecific_kcal-kg.Rmd` to analyze energy consumption of
captive spotted seals.  

    **Inputs**: `data/PHOCAS_spotted_seal_diet.csv`  
  
    **Outputs**: `data/Energy_predictions.RData`
    


8. Run integrated analysis of spotted seal consumption via `integrated_consumption_analysis.R`

    **Inputs**: `data/Combined_spotted_winter_cell_prob.tif`, `data/Combined_spotted_summer_cell_prob.tif`,
    `data/Stable_stage_matrix.RData`, `data/Seal_wt_age.RData`, `data/Energy_predictions.RData`,
    `data/prey_composition_g_kcal.RData`
    
     **Outputs**: Final results
            -
  





