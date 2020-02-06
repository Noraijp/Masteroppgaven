clear all; clc;


% make some gamma ray energies for plotiing
xp= 30:2000;

% Select efficiency calibration matrix
 effcal = '../eff_room131_10.mat'
% effcal = 'eff_room131_10.mat'
% effcal = 'eff_room131_18.mat'
%effcal = 'eff_room131_22.mat'

[efficiency, unc_efficiency, ~, ~] = efficiency_calibration(xp, effcal);


%Plotting uncertainty band
% plot(xp, efficiency,  xp, efficiency + unc_efficiency, '--', xp, efficiency - unc_efficiency, '--')
% xlabel('Gamma-Ray Energy (keV)')
% ylabel('Detector  Efficiency')
% 
% % OR

% Plotting percent uncertainty
plot(xp, 100.*unc_efficiency./ efficiency)
xlabel('Gamma-Ray Energy (keV)')
ylabel('Percent Uncertainty in Efficiency')