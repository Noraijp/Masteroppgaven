clear all; clc;

% Load parsed matrices from Zn_spec_anal.m
% load al_mat_33MeV.mat;
% load zr_mat_33MeV.mat;
% load zn_mat_33MeV.mat;
% load  y_mat_33MeV.mat
% load in_mat_33MeV.mat;
% % 
load al_mat_16MeV.mat;
load zr_mat_16MeV.mat;
load zn_mat_16MeV.mat;
load  y_mat_16MeV.mat;
load in_mat_16MeV.mat;


% Matrix format:
%
% Each row corresponds to a key energy, from the ZZZ_key_energies matrices, 
%  or from your energy assignments list
%
% Column 1:  the key energy (in data_mat)
% Column 2,10,18,..., 8.*(i-1) +2 : net gamma counts (efficiency and attenuation corrected) 
% Column 3,11,19,..., 8.*(i-1) +3 : error (in %) in net gamma counts
% Column 4,12,20,..., 8.*(i-1) +4 : sample mass #
% Column 5,13,21,..., 8.*(i-1) +5 : spectrum count time (in Julian time)
% Column 6,14,22,..., 8.*(i-1) +6 : spectrum count length (in seconds of
% live time
% Column 7,15,23,..., 8.*(i-1) +7 : parent decay half-life (in s)
% Column 8,16,24,..., 8.*(i-1) +8 : gamma branching ratio (in %)
% Column 9,17,25,..., 8.*(i-1) +9 : uncertainty in gamma b.r. (in %)

%             1XX        2XX      3XX
zr_rhodrs = mean([0.636     0.635]);
zn_rhodrs = mean([0.663     0.646]);
al_rhodrs = mean([0.197     0.196]);
in_rhodrs = mean([0.358     0.354]);
y_rhodrs =  mean([0.402     0.426]);






% Set foils for analysis
data = zn_mat;
rhodrs = zn_rhodrs;
mu_attenuation = load('./zn_xcom.txt');
pp2 = pchip(mu_attenuation(:,1).*1e3,mu_attenuation(:,2));


num_spec = (size(data,2)-1)/8;




%16MeV Zinc foils:
rows_67Cu_zn16 = [1, 2, 3];
rows_64Cu_zn16 = [11];
rows_65Ni_zn16 = [4, 10, 14];
rows_62Zn_zn16 = [6];
rows_24Na_zn16 = [12, 15];
rows_69mZn_zn16 = [5];
rows_63Zn_zn16 = [7, 9];
rows_56Mn_zn16 = [8];
rows_40K_zn16 = [13];
rows_65Zn_zn16 = [16];
data(16,:) = data(10,:);

% 33MeV Zinc foils:
% rows_67Cu_zn33 = [2, 3, 4, 7, 9];
% rows_62Zn_zn33 = [1, 5, 6, 12, 13];
% rows_65Ni_zn33 = [8, 17, 21];
% rows_69mZn_zn33 = [10];
% rows_63Zn_zn33 = [11, 14, 15, 20];
% rows_66Cu_zn33 = [16];
% rows_64Cu_zn33 = [18];
% rowa_24Na_zn33 = [19, 22];

% 16MeV zirconium foils:
% rows_97Nb_zr16 = [1];
% rows_95Nb_zr16 = [5];
% rows_95Zr_zr16 = [2, 4];
% rows_97Zr_zr16 = [3];
% rows_98Zr_zr16 = [6];
% % 
% 33MeV zirconium foils:
% rows_90mY_zr33 = [1, 4, 8];
% rows_93Y_zr33 = [2, 15];
% rows_92Y_zr33 = [3, 6, 12, 14, 17];
% rows_91mY_zr33 = [5];
% rows_97Nb_zr33 = [7];
% rows_95Zr_zr33 = [9, 11];
% rows_91Sr_zr33 = [10];
% rows_89Zr_zr33 = [13, 18, 19, 20, 21];
% rows_24Na_zr33 = [16, 22];
% 
% % 16MeV Yttrium foils:
% rows_88Y_Y16 = [1, 2, 3];
% 
% % 33MeV Yttrium foils:
% rows_90mY_y33 = [1, 3];
% rows_87mY_y33 = [2];
% rows_87Y_y33 = [4];
% rows_88Y_y33 = [5, 6]; %7
% rows_24Na_y33 = [8];

% % 16MeV Indium foils:
% rows_116mIn_in16 = [1,3, 5, 6, 9, 10, 11, 12, 13, 14];
% rows_115mIn_in16 = [2];
% rows_114mIn_in16 = [7, 8];
% rows_113mIn_in16 = [4];
% 
% % 33MeV Indium foils:
% rows_116mIn_In33 = [1, 5, 6, 7, 9, 11, 12, 14, 15, 16, 18, 19, 20, 21, 22, 24, 25, 26];
% rows_112mIn_In33 = [2];
% rows_111In_In33 = [3, 4];
% rows_115mIn_In33 = [8];
% rows_113mIn_In33 = [10];
% rows_114mIn_In33 = [13, 17];
% rows_24Na_In33 = [23, 27];
% 
% % 16MeV Aluminum foils:
% rows_24Na_al16 = [1, 2];
% 
% % 33MeV Aluminum foils:
% rows_24Na_al33 = [1, 2];



% Select rows to plot
% varToStr = @(x) inputname(1);
rows = rows_67Cu_zn16;
outName = './csv/67Cu_zn16MeV';
% rows = 12;

% Gate on foil energy - energy = 0 searches on all energies
% energy = 0;
% energy = 913;

% loop over all energies
% for energy = 0:0   % Show all foils in one plot (not for analysis!)
% for energy = 140:100:140   % Just Zirconium 16 MeV
% for energy = 240:100:240   % Just Zirconium 33 MeV
 for energy = 130:100:130   % Just Zinc 16 MeV
% for energy = 230:100:230   % Just Zinc 33 MeV
% for energy = 230:100:230   % Test mode
% for energy = 139:100:139   % Just Yttrium 16 MeV
% for energy = 239:100:239   % Just Yttrium 33 MeV
% for energy = 113:100:113   % Just Aluminum 16 MeV
% for energy = 213:100:213   % Just Aluminum 33 MeV
% for energy = 149:100:149   % Just Indium 16 MeV
% for energy = 349:100:349   % Just Indium 33 MeV


if energy==0
    masses = data(rows,8.*((1:num_spec)-1) +4);
    gammas = data(rows,8.*((1:num_spec)-1) +2);
    branching = data(rows,8.*((1:num_spec)-1) +8);
    err_gam = data(rows,8.*((1:num_spec)-1) +3);
    err_br  = data(rows,8.*((1:num_spec)-1) +9);
    
    ct = data(rows,8.*((1:num_spec)-1) +5);
    cl = data(rows,8.*((1:num_spec)-1) +6);
    lifetimes = data(rows,8.*((1:num_spec)-1) +7);
else
    [~,col] = find((data(rows,8.*((1:num_spec)-1) +4))==energy);
    % Extract unique columns
    col = unique(col');
    
%     8.*((col)-1) +4
    
    masses = data(rows,8.*(col-1) +4);
    gammas = data(rows,8.*(col-1) +2);
    branching = data(rows,8.*(col-1) +8);
    err_gam = data(rows,8.*(col-1) +3);
    err_br  = data(rows,8.*(col-1) +9);
    
    ct = data(rows,8.*(col-1) +5);
    cl = data(rows,8.*(col-1) +6);
    lifetimes = data(rows,8.*(col-1) +7);

end


lambdas = (log(2)./lifetimes);
I_gammas = branching.*.01;  %  in decimal, branching is in percent




% EoB Time
% 33 MeV
% delta_ts = (ct - juliandate(datetime('14-Aug-2018 19:46:00')) ) .* 24;
% 16 MeV
 delta_ts = (ct - juliandate(datetime('12-Aug-2018 11:00:00')) ) .* 24;

% Note: 'gammas' contains the net peak area, corrected for efficiency and
% self-attenuation
activities = gammas.*lambdas./(1.* I_gammas.*(1-exp(-lambdas.*cl) ));

% errors = sqrt(  (err_gam.* gammas./ (30 .* branching)).^2 + ...
%          ( (gammas .* err_br) ./ (1 .* (branching.^2) )  ).^2);

err_counts = lambdas.*gammas.*(.01.*err_gam)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );
err_efficiency = lambdas.*gammas.*(.04)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );
err_I_gammas = lambdas.*gammas.*(err_br./branching)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );
err_live_time = lambdas.^2.*gammas.*(5)./( (I_gammas).*4 .*( sinh(lambdas.*cl.*0.5) ).^2  );
% err_attenuation = lambdas.*gammas.*(0.05.*ppval(pp2,data(rows,1))).* (0.5 .* rhodrs)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );  
number_of_columns = size(lambdas);
number_of_columns = number_of_columns(2);
err_attenuation = lambdas.*gammas.*(0.05.*repmat(ppval(pp2,data(rows,1)),1,number_of_columns)).* (0.5 .* rhodrs)./( (I_gammas) .*(1-exp(-lambdas.*cl) )  );  

errors = sqrt( err_counts.^2 + err_efficiency.^2 + err_I_gammas.^2 + err_live_time.^2 + err_attenuation.^2);
approximate_error = activities.*sqrt((0.01.*err_gam).^2 + (0.04).^2 + (err_br./branching).^2);

% Plot each time step as same color
% errorbar(delta_ts,activities,errors,'.')
% Plot each gamma line as same color
errorbar(delta_ts',activities',errors','.')
max_time = max([max(delta_ts) 0.1]);
% xlim([0, max_time])
% Plot each time step as same color
% plot(delta_ts',activities','.')
xlabel('Time since EoB [h]')
ylabel('Activity [Bq]')
% Dump to csv for python / gnuplot


outfile = [reshape(delta_ts,[1 numel(delta_ts)])' ,reshape(activities,[1 numel(activities)])', reshape(errors,[1 numel(errors)])'];
% Turn this line on to write files out!
% csvwrite([outName '_' num2str(energy) '.dat'],outfile);

end




% write out liens of each plotted matrix into a single vector, cocatentate them into a ZZZ x 3 matreix, wrote to disk and fit in gnuplot for 
