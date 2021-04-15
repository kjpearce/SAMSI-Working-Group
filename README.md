# SAMSI-Working-Group
Codes related to global sensitivity analysis 

Sensitivity Comparison

Main file: 
rreSensComparison.m : compares parameters sensitivities computed by complex-step against sensitivity equations for enzyme cooperativity model with 2 intermed coefficients

Functions called:
rre_kinetics.m : Alen's function for RHS of enzyme cooperativity model
complex_rre_kinetics.m : RHS for complex-step 
rre_senseq.m : RHS for sensitivity equations
