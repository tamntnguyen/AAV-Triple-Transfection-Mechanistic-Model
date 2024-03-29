# AAV-Triple-Transfection-Mechanistic-Model
Mechanistic model for AAV viral vector production via triple transfection. Published June 15th 2021

This repository provides the mechanistic model in '[Mechanistic model for production of recombinant adeno-associated virus via triple transfection of HEK293 cells](https://www.cell.com/molecular-therapy-family/methods/fulltext/S2329-0501(21)00072-3)' by T.N.T Nguyen et al. The code includes the main mechanistic ode model and an example of how to simulate the triple transfection experiment. Detailed model description can be found in the article.

* **sim_full_params_mL.m** is the main file that performs the simulation of triple transfection experiment. This will result in Figure 4 in the published article

* **simE3.m** describes the triple transfection experiment setting, with media exchange at 6h post-transfection

* **ode_viralProd.m** contains the ode model for plasmid delivery and viral production* 

* **calc_growth_rate.m** computes time-dependent cell density and growth rate

* **model_params.mat**  contains optimized model parameter values



