clear all; clc;

% Fitted efficiency parameters from gf3
a_22cm = -0.192889;
b_22cm = 1.46672;
c_22cm = -0.0256788;

% efficiency function has form 
% exp(a*log(x)**2+b*log(x)+c)

eff_22cm = @(x)(exp(a_22cm.*(log(x)).^2+b_22cm.*log(x)+c_22cm));

energy = 898.042;

eff_22 = eff_22cm(energy)



a_10cm = -0.206303;
b_10cm = 1.6182;
c_10cm = -1.18783;

% efficiency function has form 
% exp(a*log(x)**2+b*log(x)+c)

eff_10cm = @(x)(exp(a_10cm.*(log(x)).^2+b_10cm.*log(x)+c_10cm));

energy = 898.042;

eff_10 = eff_10cm(energy)


% % Calculate flux
    
N_gammas_monitor = [594 20028];                             % gammas
lambda_monitor = log(2)./[99.476*60 4.486*3600];            % 1/s 
Nt_monitors = 0.550.*6.022E23.*[0.0428 0.9572]./114.818;    % nuclei
eff_monitors = [eff_10cm(391.698) eff_10cm(336.241)];       % percent
I_gamma_monitors = [0.6494 0.459];                          % percent
t_i = (8*3600) - (3*60);                                    % s
t_d = (2*3600) + (54*60) +36;                               % s
t_c = 201;                                                  % s
sigma_monitors = [200 200].*1E-27;                          % cm^2

flux_monitors = (N_gammas_monitor.* lambda_monitor)./(sigma_monitors.*Nt_monitors.*eff_monitors.*I_gamma_monitors.*(1-exp(-lambda_monitor.*t_i)).*exp(-lambda_monitor.*t_d).*(1-exp(-lambda_monitor.*t_c)))

N_gammas_target = [1039 1007];                             % gammas
lambda_target = log(2)./[12.701*3600 61.83*3600];          % 1/s 
Nt_targets = 0.8463333.*6.022E23./65.38;                   % nuclei
eff_targets = [eff_10cm(1345.77) eff_10cm(184.577)];       % percent
I_gamma_targets = [0.00475 0.487];                         % percent
t_i = (8*3600) - (3*60);                                   % s
t_d = (1*3600) + (19*60) +56;                              % s
t_c = 5510;                                                % s
sigma_targets = [200 200].*1E-27;                          % cm^2

xs_targets = (N_gammas_target.* lambda_target)./(flux_monitors.*Nt_targets.*eff_targets.*I_gamma_targets.*(1-exp(-lambda_target.*t_i)).*exp(-lambda_target.*t_d).*(1-exp(-lambda_target.*t_c)))./1E-27
