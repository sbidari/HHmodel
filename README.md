# Codes to run the model in 'The Impact of Household Sizes on Measles Transmission: A Long-term Perspective'
ClassMixing.m and dependencies (State_to_Class_Mat_With_Elders.m,State_to_Class.m) - Generates mixing data for the models using data from dataInput folder for each year<br>
contactData.m - create contact times data used by the model using POLYMOD data<br>

Get_r0_seir.m - calculates R0 using early growth rate<br>
convert_tau_to_R0.m - convert unit time transmission rate to basic reproduction number (R0)<br>

run_epi_1968_2019.m - Generate yearly infections from 1950-2000<br>
run_epi_1968_2019_fixed_hh.m - Generate yearly infections from assuming no HH size changes from 1970. <br>
These are called by run_epi_*<br>
Get_Q_seir.m - Generate epidemic transition matrix<br>
Get_Qdemo_seir.m - Generate demographic transition matrix<br>
Get_Eq_Demography_seir.m - Generate Equilibrium HH states<br>
Get_I_seir_deterministic/stochastic.m - Solve the ode for new infections<br>
Get_initial_inf_seir.m - Obtain initial conditions for the ode model<br>
create_brand.cpp - mex code for genrating binomial random variable <br>
