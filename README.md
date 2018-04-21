# Joint Modelling for Organ Transplantation Outcomes for Patients with Diabetes and the End-Stage Renal Disease

We have included the simulation codes for three simulation studies in three folders, respectively. 

For each subfolder,  

Source('main.R') to run simulations.

For the 9 .R files in the folder:

1. 'Simulation.R' is used to generating simulation data. If we want to try different simulation settings, we can modify this file.

2. 'Estep.R' contains functions for the E-step of MCEM step. 

3. 'Mstep.R' contains functions for the M-step of MCEM step. 

4. 'InitialValue.R' contains functions for obtaining good initial value for MCEM. We fit the mixed effects model and survival model seperated to obtain good initial values.

5. 'optimizing.R' contains function for performing E-step and M-step iteratively.

6. 'main.R' is used to run the whole algorithm. Set 'DoBootstrap = 1' if you want to perform bootstrap method. 

7. 'Bootstrap.R' contains functions for performing bootstrap method. 

8. './bootstrap/Initial_bootstrap.R' contains functions for obtaining good initial value for MCEM in bootstrap.

9. './bootstrap/optimization_bootstrap.R' contains functions for E-step and M-step in bootstrap.

