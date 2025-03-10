function [Z_real, S11_calc, S21_calc] =  ...
    input_imped_calc(freq, abs_pwr_ant_red, abs_pwr_bckd, Z0_wire, gamma_diel_slab, t_slab, eps_diel)
%% This function calculates the input impedance from the extinction power using
% the equation (9). Then, assuming the two-port network model of the EELS
% experiment, it calculates S-parameters using its input impedance throught
% the more convenient ABCD matrices multiplication
%% Parameters:
%   freq - frequnency mesh
%   abs_pwr_ant_red - extinction power reduced by those calculated for the
%                     background (length should be the same as for "freq" 
%                     variable)
%   Z0_wire - characteristic impedance of the electron beam, theoretically
%             evaluated as a TM-wave impedance with k_el = omega/v_el, or 
%             directly simulated (our option)
%   gamma_diel_slab - propagation constant (alpha + i*beta) of the electron
%                     beam propagating in the silicon nitride dielectric 
%                     membrane. We simulated it and imported as a separate
%                     file
%   t_slab - silicon nitride membrane thickness (100 nm, as an example)
%   eps_diel - complex dielectric constant (eps_real + i*eps_imag) of the
%              silicon nitride membrane (taken from reference [70])

phys_const_generate

omega = 2.*pi.*freq;
k0 = omega.*sqrt(eps_0.*mu_0);

Z0_free = sqrt(mu_0./eps_0).*ones(size(freq));                                          % Free space input impedance (120*pi)
Z0_slab = Z0_free./sqrt(eps_diel).*(-1j).*gamma_diel_slab./(k0.*sqrt(eps_diel));        % Characteristic impedance of the electron beam in the silicon nitride membrane 
Y0_slab = 1./Z0_slab;                                                                   % Slab admittance

%% Input impedance calculation
% ABCD parameters for the transmission line represented by the electron
% beam traversing the silicon nitride membrane
A = cosh(gamma_diel_slab.*t_slab);
B = Z0_slab.*sinh(gamma_diel_slab.*t_slab);
C = Y0_slab.*sinh(gamma_diel_slab.*t_slab);
D = cosh(gamma_diel_slab.*t_slab);
delta_S = A + B./Z0_wire + C.*Z0_wire + D;
xi = A + B./Z0_wire;

% Equation (9)
alpha = ((abs_pwr_ant_red + abs_pwr_bckd).*(abs(delta_S).^2.*abs_pwr_bckd + 4) - ...
    abs(delta_S).^2.*abs_pwr_bckd)./(4.*abs(xi).^2.*(1 - abs_pwr_ant_red - abs_pwr_bckd));
Z_real = Z0_wire./alpha;

%% S-parameters calculation
% ABCD parameters for the antenna represented as the admittance connected
% in parallel to the electron beam
Y_real = 1./Z_real;
A_ant = ones(size(freq));
B_ant = zeros(size(freq));
C_ant = Y_real;
D_ant = ones(size(freq));

for ind = 1:length(freq)
    A_ant_mtx = [A_ant(ind), B_ant(ind); C_ant(ind), D_ant(ind)];
    A_sub_mtx = [A(ind), B(ind); C(ind), D(ind)];
    A_full = A_ant_mtx*A_sub_mtx;

    A_fin(ind, 1) = A_full(1,1);
    B_fin(ind, 1) = A_full(1,2);
    C_fin(ind, 1) = A_full(2,1);
    D_fin(ind, 1) = A_full(2,2);

end

% Conversion from ABCD parameters to S-parameters for the entire two-port network 
[S11_calc, S12_calc, S21_calc, S22_calc] = ABCDtoSparams(A_fin, B_fin, C_fin, D_fin, Z0_wire);

end