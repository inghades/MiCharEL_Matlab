
%% Pure background-substrate only S-parameters (COMSOL simulated)
fileName_S11_bckd_real = "data/S11ant_bckd_real.txt";
fileName_S11_bckd_imag = "data/S11ant_bckd_imag.txt";
fileName_S21_bckd_real = "data/S21ant_bckd_real.txt";
fileName_S21_bckd_imag = "data/S21ant_bckd_imag.txt";
%% Wire propagation constant and characteristic impedances (simulated)
fileName_Z0 = "data/Z0_const_wire_real.txt";
fileName_beta = "data/beta_prop_const_wire.txt";                                 % Beam in vacuum (wavenumber)
fileName_alpha_wire_diel = "data/alpha_prop_const_wire_diel.txt";                % Beam in dielectric (absorption)
fileName_beta_wire_diel = "data/beta_prop_const_wire_diel.txt";                  % Beam in dielectric (wavenumber)
%% Antenna propagation constant (simulated)
fileName_alpha_ant = "data/alpha_prop_const_ant_3D.txt";                         % Plasmonic dipole (strip) mode (absorption)
fileName_beta_ant = "data/beta_prop_const_ant_3D.txt";                           % Plasmonic dipole (strip) mode (wavenumber)
fileName_eps = "data/eps_real_SN40_300.txt";                                     % Dielectric constant of the SiN membrane (real) taken from [70]
fileName_eps_imag = "data/eps_imag_SN40_300.txt";                                % Dielectric constant of the SiN membrane (imag) taken from [70]


if spat_pt_name == "edge"                                                    % Edge of the dipole

    %% Experimental spectra
    filename_real_spec = "data/EELS_probability_edge_av.txt";                    % Mean value averaged over 16x16 pixels
    filename_real_spec_plus = "data/EELS_probability_edge_av_plus.txt";          % Mean value plus standard deviation (std)
    filename_real_spec_minus = "data/EELS_probability_edge_av_minus.txt";        % Mean value minus standard deviation (std)
    %% Antenna S-parameters (COMSOL Simulated)
    fileName_S11_ant_real = "data/S11ant_edge_real_dist40nm.txt";
    fileName_S11_ant_imag = "data/S11ant_edge_imag_dist40nm.txt";
    fileName_S21_ant_real = "data/S21ant_edge_real_dist40nm.txt";
    fileName_S21_ant_imag = "data/S21ant_edge_imag_dist40nm.txt";
    %% Antenna input impedance (COMSOL simulated)
    filename_Zin_real = "data/Zin_ant_edge_re.txt";
    filename_Zin_imag = "data/Zin_ant_edge_im.txt";
    
    ind_res = [26; 94; 158; 219];                                           % Indices for resonances corresponding to frequency mesh

elseif spat_pt_name == "interm"                                             % Intermediate point (corresponding to L_dipole/3 from the edge)

    %% Experimental spectra
    filename_real_spec = "data/EELS_probability_two_third_av.txt";               % Mean value averaged over 16x16 pixels
    filename_real_spec_plus = "data/EELS_probability_two_third_av_plus.txt";     % Mean value plus standard deviation (std)
    filename_real_spec_minus = "data/EELS_probability_two_third_av_minus.txt";   % Mean value minus standard deviation (std)
    %% Antenna S-parameters (COMSOL Simulated)
    fileName_S11_ant_real = "data/S11ant_interm_real_dist40nm.txt";
    fileName_S11_ant_imag = "data/S11ant_interm_imag_dist40nm.txt";
    fileName_S21_ant_real = "data/S21ant_interm_real_dist40nm.txt";
    fileName_S21_ant_imag = "data/S21ant_interm_imag_dist40nm.txt";
    %% Antenna input impedance (COMSOL simulated)
    filename_Zin_real = "data/Zin_ant_interm_re.txt";
    filename_Zin_imag = "data/Zin_ant_interm_im.txt";
    
    ind_res = [158];                                                        % Indices for resonances corresponding to frequency mesh

elseif spat_pt_name == "center"

    %% Experimental spectra
    filename_real_spec = "data/EELS_probability_center_av.txt";                  % Mean value averaged over 16x16 pixels
    filename_real_spec_plus = "data/EELS_probability_center_av_plus.txt";        % Mean value plus standard deviation (std)
    filename_real_spec_minus = "data/EELS_probability_center_av_minus.txt";      % Mean value minus standard deviation (std)
    %% Antenna S-parameters (COMSOL Simulated)
    fileName_S11_ant_real = "data/S11ant_center_real_dist40nm.txt";
    fileName_S11_ant_imag = "data/S11ant_center_imag_dist40nm.txt";
    fileName_S21_ant_real = "data/S21ant_center_real_dist40nm.txt";
    fileName_S21_ant_imag = "data/S21ant_center_imag_dist40nm.txt";
    %% Antenna input impedance (COMSOL simulated)
    filename_Zin_real = "data/Zin_ant_center_re.txt";
    filename_Zin_imag = "data/Zin_ant_center_im.txt";

    ind_res = [94; 219];                                                    % Indices for resonances corresponding to frequency mesh

end

%% Simulatied characteristics import
% S-parameters for a background
[freq, S11_bckd_real] = importfile(fileName_S11_bckd_real);
[~, S11_bckd_imag] = importfile(fileName_S11_bckd_imag);
[freq_bckd, S21_bckd_real] = importfile(fileName_S21_bckd_real);
S21_bckd_real = spline(freq_bckd, S21_bckd_real, freq);
[freq_bckd, S21_bckd_imag] = importfile(fileName_S21_bckd_imag);
S21_bckd_imag = spline(freq_bckd, S21_bckd_imag, freq);

% S-parameters for the whole system (antenna + background substrate)
[freq_ant, S11_ant_real] = importfile(fileName_S11_ant_real);
S11_ant_real = spline(freq_ant, S11_ant_real, freq);
[freq_ant, S11_ant_imag] = importfile(fileName_S11_ant_imag);
S11_ant_imag = spline(freq_ant, S11_ant_imag, freq);
[freq_ant, S21_ant_real] = importfile(fileName_S21_ant_real);
S21_ant_real = spline(freq_ant, S21_ant_real, freq);
[freq_ant, S21_ant_imag] = importfile(fileName_S21_ant_imag);
S21_ant_imag = spline(freq_ant, S21_ant_imag, freq);

% Wire microwave parameters (Z0, beta)
[freq_Z0, Z0_wire] = importfile(fileName_Z0);
Z0_wire = spline(freq_Z0, Z0_wire, freq);
[freq_beta, beta_wire] = importfile(fileName_beta);
beta_wire = spline(freq_beta, beta_wire, freq);

% Wire in dielectric microwave parameters (Z0, beta)
[freq_beam, alpha_wire_diel] = importfile(fileName_alpha_wire_diel);
alpha_wire_diel = spline(freq_beam, alpha_wire_diel, freq);
[freq_beam, beta_wire_diel] = importfile(fileName_beta_wire_diel);
beta_wire_diel = spline(freq_beam, beta_wire_diel, freq);
gamma_wire_diel = alpha_wire_diel + 1j.*beta_wire_diel;

%% Antenna microwave parameters (alpha, beta)
[freq_ant, alpha_ant] = importfile(fileName_alpha_ant);
alpha_ant = spline(freq_ant, alpha_ant, freq);
[freq_ant, beta_ant] = importfile(fileName_beta_ant);
beta_ant = spline(freq_ant, beta_ant, freq);
gamma_ant = alpha_ant + 1j.*beta_ant;

% Input impedance (COMSOL simulated)
[freq_Z, Ztrue_in_re] = importfile(filename_Zin_real);                     % This defines the frequency range
Ztrue_in_re = spline(freq_Z, Ztrue_in_re, freq);
[freq_Z, Ztrue_in_im] = importfile(filename_Zin_imag);                     % This defines the frequency range
Ztrue_in_im = spline(freq_Z, Ztrue_in_im, freq);

% Substrate parameters (eps_r, eps_i)
[freq_eps, eps_diel] = importfile(fileName_eps);
eps_diel = spline(freq_eps, eps_diel, freq);
[freq_eps_i, eps_diel_extend_i] = importfile(fileName_eps_imag);
eps_diel_i = spline(freq_eps_i, eps_diel_extend_i, freq);
eps_diel = eps_diel - 1j.*eps_diel_i;

S11_ant = S11_ant_real + 1j.*S11_ant_imag;
S21_ant = S21_ant_real + 1j.*S21_ant_imag;
S11_bckd = S11_bckd_real + 1j.*S11_bckd_imag;
S21_bckd = S21_bckd_real + 1j.*S21_bckd_imag;

S11sqr_ant = abs(S11_ant).^2;
S21sqr_ant = abs(S21_ant).^2;
S11sqr_bckd = abs(S11_bckd).^2;
S21sqr_bckd = abs(S21_bckd).^2;

%% Experimental spectra
[freq_vac, eels_spec_exper] = importfile(filename_real_spec);                        % Experimental EEL spectrum at one single point along the dipole
Spec_eels_exper = spline(freq_vac, eels_spec_exper, freq);
[freq_plus, eels_spec_exper_plus] = importfile(filename_real_spec_plus);                  % Experimental EEL spectrum at one single point along the dipole (+ sigma)
Spec_eels_exper_plus = spline(freq_plus, eels_spec_exper_plus, freq);
[freq_minus, eels_spec_exper_minus] = importfile(filename_real_spec_minus);                 % Experimental EEL spectrum at one single point along the dipole (- sigma)
Spec_eels_exper_minus = spline(freq_minus, eels_spec_exper_minus, freq);