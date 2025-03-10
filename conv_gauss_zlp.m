function [conv_spec] = conv_gauss_zlp(freq, en_fwhm, init_eel_spec)

qe = 1.6e-19;
h_bar = 1.055e-34;
h_bareV = h_bar./qe;

fwhw_freq = (en_fwhm./2)./h_bareV./(2.*pi);
sigma = fwhw_freq./sqrt(2.*log(2));         % Dispersion of the Gaussian Kernel for FWHM

delt_freq = freq(2) - freq(1);
delt_energy = delt_freq.*h_bareV.*2.*pi;
freq_ZLP = (0:delt_freq:(delt_freq.*(2.*length(freq) - 1))).' - (delt_freq.*(2.*length(freq) - 1))./2;

%% Gaussian Kernel
gauss_kernel = 1./sqrt(2.*pi.*sigma.^2).*exp(-1./2.*(freq_ZLP./sigma).^2);
bckd_one_arm_ZLP_new = gauss_kernel.*delt_freq./delt_energy;
conv_spec = conv(init_eel_spec, bckd_one_arm_ZLP_new.*delt_energy, 'same');

end