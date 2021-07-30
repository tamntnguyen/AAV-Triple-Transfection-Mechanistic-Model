# AAV-Triple-Transfection-Mechanistic-Model
Mechanistic model for AAV viral vector production via triple transfection. Published June 15th 2021


The code in this repository includes the ode model used in 'Mechanistic model for production of recombinant adeno-associated virus via triple transfection of HEK293 cells' by T.N.T Nguyen et al.

sim_full_params_mL.m - Main file to execute the triple transfection experiment. This will result in Figure 4 in the published article
simE3.m              - File that describes the triple transfection experiment setting, with media exchange at 6h post-transfection
ode_viralProd.m      - File that contains the ode model for plasmid delivery and viral production
calc_growth_rate.m   - Funciton to compute time-dependent cell density and growth rate
model_params.mat     - Optimized model parameters



