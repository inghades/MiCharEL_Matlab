%% Short description of the code
The code is used to convert the EEL probability to the input impedance and S-parameters of the plasmonic antenna (or any other plasmonic structure). The code uses the algorithm illustrated in Fig. 4 of the manuscript.

To run the script:
-- open the main script file named main_convolved_spectrum_conversion.m in MATLAB software (code was developed using R2022b version)
-- run the script

By default, the script uses the experimental spectrum collected at the center of the plasmonic dipole and gives the extracted microwave parameters as an output. The identical results of the algorithm performance is shown in the right column of Fig.5 in the main manuscript.

To run the algorithm for the plasmonic dipole edge:
-- delete 'center' in the spat_pt_name variable at line 10 of the main script
-- instead type 'edge' in the spat_pt_name variable
-- run the main script


