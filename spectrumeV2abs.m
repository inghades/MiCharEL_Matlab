function [pwr_abs_cnvtd] = spectrumeV2abs(freq, eels_spectrum_eV, r_beam, E_beam)
%% This spectrum converts the EEL spectrum into the extinction power by 
% applying the reverse function to equation (4)

phys_const_generate

vel = sqrt(2.*E_beam.*qe./me);
gammaL = 1./sqrt(1 - vel.^2./c0.^2);

omega = 2.*pi.*freq;
var_bessel = omega.*r_beam./(vel.*gammaL);

area_integr_part = r_beam.^2.*(besselk(1, var_bessel).^2 - besselk(0, var_bessel).*besselk(2, var_bessel));
power_const = -2.*pi.*(k.^2.*qe.^2.*omega.^2)./(vel.^3.*gammaL.^2).*eps_0;

pwr_abs_cnvtd = eels_spectrum_eV.*h_bareV.*(pi.*h_bar.*omega)./2;                     % Reverse function to the equation (4)
pwr_abs_cnvtd = pwr_abs_cnvtd./(power_const.*area_integr_part);

end

