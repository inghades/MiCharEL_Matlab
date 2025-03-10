%% The script is created by Igor Getmanov (e-mail igor.getmanov@kaust.edu.sa)
%% Version 1.0

clear all, close all

phys_const_generate

%% User Configurable parameters (need different simulated data for another beam energy & radius)
E_beam = 80e3;                                                              % Beam energy
r_beam = 2.5e-9;                                                            % Beam radius
bin_factor = 16;                                                            % Energy binning
t_sub = 100e-9;                                                             % Silicon Nitride substrate thickness
spat_pt_name = 'center';                                                    % Name of the point presented in paper
coeff = 1.58;                                                               % Correction coefficient for the ZLP effect (Supp. Info V)
%% Supplimentary constants
v = sqrt(2.*E_beam.*qe./me);                                                % Electron velocity
gamma_L = 1./sqrt(1 - v.^2./c0.^2);                                         % Lorentz factor

%% Gaussian Kernel parameters (ZLP kernel for convolution of simulated spectra)
mu = 0;                                                                     % Mean value
x_fwhm = 0.06;                                                              % FWHM of the ZLP

import_data                                                                 % Import data (background parameters & other simulated parameters for comparison)     

omega = 2.*pi.*freq;                                                        % Angular frequency mesh based on improted data from import_data
k0 = omega.*sqrt(eps_0.*mu_0);                                              % Free space wavenumber

%% Extinction power for antenna + background and background alone calculated using corrected approach
abs_pwr_ant = (1 - S11sqr_ant - S21sqr_ant)./(1 - S11sqr_ant);              % Equation (8) for antenna + substrate
abs_pwr_bckd = (1 - S11sqr_bckd - S21sqr_bckd)./(1 - S11sqr_bckd);          % Equation (8) for substrate alone

figure(1),
plot(freq./1e12, abs_pwr_ant, 'r', 'DisplayName', "Antenna + Substrate")
hold on
plot(freq./1e12, abs_pwr_bckd, 'k', 'DisplayName', "Substrate alone")
title("Extinction power. Corrected approach")
xlabel("Frequency (THz)")
ylabel("P_{ext}/P_0")
xlim([25,150])
legend
grid on
saveas(gcf,['Extinction_power_Corrected_approach_', spat_pt_name, '.png'])

%% Convertion of the simulated extinction power to EEL spectra (equation (4))
[~, eels_spect_eV_sim] = abs2spectrum(freq, abs_pwr_ant, r_beam, E_beam);            % Simulated spectrum (Antenna + substrate)
[~, eels_spect_eV_bckd_sim] = abs2spectrum(freq, abs_pwr_bckd, r_beam, E_beam);      % Simulated spectrum (Substrate alone)

figure(2),
plot(freq./1e12, eels_spect_eV_sim, 'r', 'DisplayName', "Antenna + Substrate")
hold on
plot(freq./1e12, eels_spect_eV_bckd_sim, 'k', 'DisplayName', "Substrate alone")
plot(freq./1e12, eels_spect_eV_sim - eels_spect_eV_bckd_sim, 'b', 'DisplayName', "Background reduced")
title("Simulated spectrum: pure (no ZLP spreading)")
xlabel("Frequency (THz)")
ylabel("EEL probability (eV^{-1})")
xlim([25,150])
legend
grid on
saveas(gcf,['Simulated_spectrum_no_ZLP_spread_', spat_pt_name, '.png'])

%% Convolution with the gaussian ZLP of 60 meV energy resolution
[Spec_conv_eels_gauss] = conv_gauss_zlp(freq, x_fwhm, eels_spect_eV_sim);            % Convoluion of the simulated spectrum with the gaussian ZLP (antenna + substrate)
[Spec_conv_eels_bckd_gauss] = conv_gauss_zlp(freq, x_fwhm, eels_spect_eV_bckd_sim);  % Convoluion of the simulated spectrum with the gaussian ZLP (substrate alone)
Spec_conv_reduced = Spec_conv_eels_gauss - Spec_conv_eels_bckd_gauss;                % Reduced antenna + substracte spectrum by the background

figure(3),
plot(freq./1e12, Spec_eels_exper, 'k', 'DisplayName', "Experimental")
hold on
plot(freq./1e12, Spec_conv_reduced, 'r', 'DisplayName', "Simulated (antenna + substrate & background reduced)")
plot(freq./1e12, Spec_eels_exper_plus, 'b--', 'DisplayName', 'Deviation caused by averaging over 16 x 16 pixels')
plot(freq./1e12, Spec_eels_exper_minus, 'b--', 'HandleVisibility', 'off')
title("Simulated spectrum: ZLP spreaded effect included")
xlabel("Frequency (THz)")
ylabel("EEL probability (eV^{-1})")
xlim([25, 150])
ylim([0, 0.06])
legend
grid on
saveas(gcf,['Simulated_spectrum_ZLP_spread_effect_included_', spat_pt_name, '.png'])

%% Manually defined resonance frequneices
% Simulated spectra
[abs_pwr_ant_conv_reverse] = spectrumeV2abs(freq, Spec_conv_reduced, r_beam, E_beam);
[abs_pwr_bckd_conv_reverse] = spectrumeV2abs(freq, Spec_conv_eels_bckd_gauss, r_beam, E_beam);
% Experimental spectra
[abs_pwr_ant_reverse_exper] = spectrumeV2abs(freq, Spec_eels_exper, r_beam, E_beam);                % Experiment (averaged over 16 x 16)
[abs_pwr_ant_exper_plus] = spectrumeV2abs(freq, Spec_eels_exper_plus, r_beam, E_beam);              % Experiment (averaged + std)
[abs_pwr_ant_exper_minus] = spectrumeV2abs(freq, Spec_eels_exper_minus, r_beam, E_beam);            % Experiment (averaged - std)

figure(4),
plot(freq./1e12, abs_pwr_ant, 'k', 'DisplayName', 'Unmodified simulated extinction power')
hold on
plot(freq./1e12, abs_pwr_ant_conv_reverse, 'r', 'DisplayName', 'Extracted extinction power from simulated spectrum')
plot(freq./1e12, abs_pwr_ant_reverse_exper, 'b', 'DisplayName', 'Extracted extinction power from experimental spectrum')
title("Extracted extinction power: comparison")
xlabel("Frequency (THz)")
ylabel("Extinction power P_{ext}/P_{0}")
xlim([25, 150])
ylim([0, 0.6])
legend
grid on
saveas(gcf,['Extracted_extinction_power_comparison_', spat_pt_name, '.png'])

%% Extraction of input impedance and S-parameters from the extinction power (equation (9))
[Z_real_sim, S11_calc_sim, S21_calc_sim] =  ...                                                 % Extraction with our algorithm from the simulated spectrum
    input_imped_calc(freq, abs_pwr_ant_conv_reverse.*coeff, ...
    abs_pwr_bckd, Z0_wire, gamma_wire_diel, t_sub, eps_diel);
[Z_real_exper, S11_calc_exper, S21_calc_exper] =  ...                                           % Extraction with our algorithm from the experimental spectrum
    input_imped_calc(freq, abs_pwr_ant_reverse_exper.*coeff, abs_pwr_bckd, ...
    Z0_wire, gamma_wire_diel, t_sub, eps_diel);
% Calculations for averaging errors
[Z_real_exper_plus, S11_calc_exper_plus, S21_calc_exper_plus] =  ...                            % Extraction with our algorithm from the experimental spectrum (+ std)
    input_imped_calc(freq, abs_pwr_ant_exper_plus.*coeff, abs_pwr_bckd, ...
    Z0_wire, gamma_wire_diel, t_sub, eps_diel);
[Z_real_exper_minus, S11_calc_exper_minus, S21_calc_exper_minus] =  ...                         % Extraction with our algorithm from the experimental spectrum (- std)
    input_imped_calc(freq, abs_pwr_ant_exper_minus.*coeff, abs_pwr_bckd, ...
    Z0_wire, gamma_wire_diel, t_sub, eps_diel);
% For errorbars calculation
dZ_exp_pos = Z_real_exper_minus - Z_real_exper;
dZ_exp_neg =  Z_real_exper - Z_real_exper_plus;
dS11_exp_pos = pow2db(abs(S11_calc_exper_plus).^2) - pow2db(abs(S11_calc_exper).^2);
dS11_exp_neg = pow2db(abs(S11_calc_exper).^2) - pow2db(abs(S11_calc_exper_minus).^2);
dS21_exp_pos = pow2db(abs(S21_calc_exper_minus).^2) - pow2db(abs(S21_calc_exper).^2);
dS21_exp_neg = pow2db(abs(S21_calc_exper).^2) - pow2db(abs(S21_calc_exper_plus).^2);

figure(5),
%% Zin plot
subplot(3,1,1)
plot(freq./1e12, Ztrue_in_re./Z0_wire, 'k', 'DisplayName', "Standard Full-Wave Simulation")                                       % Z_true_in is the input impedance calculated based on COMSOL full wave simulations
hold on
plot(freq./1e12, Z_real_sim./Z0_wire, 'r--', 'DisplayName', "Extracted through algorthm using simulated spectrum")
plot(freq./1e12, abs(Z_real_exper./Z0_wire), 'r-', 'DisplayName', "Extracted through algorthm using simulated spectrum")
plot(freq(ind_res)./1e12, Ztrue_in_re(ind_res)./Z0_wire(ind_res), 'ko', 'HandleVisibility', 'off')
plot(freq(ind_res)./1e12, Z_real_sim(ind_res)./Z0_wire(ind_res), 'ro', 'HandleVisibility', 'off')
plot(freq(ind_res)./1e12, Z_real_exper(ind_res)./Z0_wire(ind_res), 'rd', 'HandleVisibility', 'off')
errorbar(freq(ind_res)./1e12, Z_real_exper(ind_res)./Z0_wire(ind_res), ...
    dZ_exp_neg(ind_res)./Z0_wire(ind_res), dZ_exp_pos(ind_res)./Z0_wire(ind_res), ...
    'bd', 'HandleVisibility', 'off', 'LineWidth', 2)
xlim([25, 150])
ylim([0, 7])
ylabel("Z_{in}/Z^{beam}_0")
legend
grid on

%% S11 plot
subplot(3,1,2)
plot(freq./1e12, pow2db(S11sqr_ant), 'k')                                                                                              % Simulated S11 through COMSOL
hold on
plot(freq./1e12, pow2db(abs(S11_calc_sim).^2), 'r--')
plot(freq./1e12, pow2db(abs(S11_calc_exper).^2), 'r-')
plot(freq(ind_res)./1e12, pow2db(S11sqr_ant(ind_res)), 'ko')
plot(freq(ind_res)./1e12, pow2db(abs(S11_calc_sim(ind_res)).^2), 'ro')
plot(freq(ind_res)./1e12, pow2db(abs(S11_calc_exper(ind_res)).^2), 'rd')
errorbar(freq(ind_res)./1e12, pow2db(abs(S11_calc_exper(ind_res)).^2), ...
    dS11_exp_neg(ind_res), dS11_exp_pos(ind_res), 'bd', 'HandleVisibility', 'off', ...
    'LineWidth', 2)
xlim([25, 150])
ylim([-18, 0])
ylabel("S_{11} (dB)")
grid on

subplot(3,1,3)
plot(freq./1e12, pow2db(S21sqr_ant), 'k')
hold on
plot(freq./1e12, pow2db(abs(S21_calc_sim).^2), 'r--')
plot(freq./1e12, pow2db(abs(S21_calc_exper).^2), 'r')
plot(freq(ind_res)./1e12, pow2db(S21sqr_ant(ind_res)), 'ko')
plot(freq(ind_res)./1e12, pow2db(abs(S21_calc_sim(ind_res)).^2), 'ro')
plot(freq(ind_res)./1e12, pow2db(abs(S21_calc_exper(ind_res)).^2), 'rd')
errorbar(freq(ind_res)./1e12, pow2db(abs(S21_calc_exper(ind_res)).^2), ...
    dS21_exp_neg(ind_res), dS21_exp_pos(ind_res), 'bd', 'HandleVisibility', 'off', ...
    'LineWidth', 2)
xlim([25, 150])
ylim([-7, 0])
ylabel("S_{21} (dB)")
xlabel("Frequency (THz)")
grid on
saveas(gcf,['Extracted_microwave_parameters_', spat_pt_name, '.png'])