function [eels_spectrum, eels_spectrum_eV] = abs2spectrum(freq, pwr_abs, r_beam, E_beam)
%% This function converts the extinction power of EM wave to the EEL spectrum in accordance with equation (6)

phys_const_generate

vel = sqrt(2.*E_beam.*qe./me);
gammaL = 1./sqrt(1 - vel.^2./c0.^2);

omega = 2.*pi.*freq;
var_bessel = omega.*r_beam./(vel.*gammaL);

area_integr_part = r_beam.^2.*(besselk(1, var_bessel).^2 - besselk(0, var_bessel).*besselk(2, var_bessel));
power_const = -2.*pi.*(k.^2.*qe.^2.*omega.^2)./(vel.^3.*gammaL.^2).*eps_0;

eels_spectrum = 2./(pi.*h_bar.*omega).*power_const.*area_integr_part.*pwr_abs;  % Equation (4)
eels_spectrum_eV = eels_spectrum./h_bareV;                                      % Convertion to dimensions (eV^-1)

end

