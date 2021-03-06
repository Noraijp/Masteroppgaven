clear all; clc;
fclose('all');


% loading key_energies and glines.
load al_key_energies.mat
load zn_key_energies.mat
load zr_key_energies.mat
load y_key_energies.mat
load in_key_energies.mat
%  
% % Gamma line half-lives and intensities
load al_glines.mat
load zn_glines.mat
load zr_glines.mat
load in_glines.mat
load y_glines.mat


% Load file names
zn_fn = textscan(fopen('zn_fnames.txt'),'%s');   
zr_fn = textscan(fopen('zr_fnames.txt'),'%s');   
y_fn = textscan(fopen('y_fnames.txt'),'%s');   
al_fn = textscan(fopen('al_fnames.txt'),'%s');   
in_fn = textscan(fopen('in_fnames.txt'),'%s');   


% Set areal densities and their uncertainties
%            1XX        2XX      3XX
zr_rhodrs = [0.636     0.635];
zn_rhodrs = [0.663     0.646];
al_rhodrs = [0.197     0.196];
in_rhodrs = [0.358     0.000     0.354];
y_rhodrs =  [0.402     0.426];

unc_zn_rhodrs = [0.024   0.011];
unc_zr_rhodrs = [0.0045  0.0036];
unc_y_rhodrs  = [0.0017  0.00053];
unc_in_rhodrs = [0.0071  0.000   0.0040];
unc_al_rhodrs = [0.00055 0.00160];

%%
clc; 

% Choose files for analysis
fitzpeaks_reports = al_fn;
key_energies = al_key_energies;
glines = al_glines;
% EoB Time
% 16 MeV
EoB_Time = '12-Aug-2018 11:00:00';
% 33 MeV
%EoB_Time = '14-Aug-2018 19:46:00';
rhodrs = al_rhodrs;
mu_attenuation = load('../al_xcom.txt');
unc_rhodrs = unc_al_rhodrs;


% Test new fitzpeaks_parser wrapper function

% data_mat  = fitzpeaks_parser(filename_list, key_energies, gamma_lines, EoB_Time, attenuation_data, rhodrs, unc_rhodrs)

%calls the fitzpeaks parser
output = fitzpeaks_parser(string(fitzpeaks_reports{1}), key_energies, glines, EoB_Time, mu_attenuation, rhodrs, unc_rhodrs);



%%

% Load parsed matrices from output of previous section
% load ni_mat.mat;
% load fe_mat.mat;
% load cu_mat.mat;
% load ir_mat.mat;


% Matrix format:
%
% Each row corresponds to a key energy, from the ZZZ_key_energies matrices,
%  or from your energy assignments list
%
% Each column (one spectrum) corresponds to the GammaCounts objects containing
% info about the peaks in that spectrum



% Set foils for analysis
data = output;


% Rows for each decay product:
%16MeV Zinc foils:
% rows_67Cu_zn16 = [1, 2, 3];
% rows_64Cu_zn16 = [11];
% rows_65Ni_zn16 = [4, 10, 14];
% rows_62Zn_zn16 = [6];
% rows_24Na_zn16 = [12, 15];
% rows_69mZn_zn16 = [5];
% rows_63Zn_zn16 = [7, 9];
% rows_56Mn_zn16 = [8];
% rows_40K_zn16 = [13];
% rows_65Zn_zn16 = [16]; %sp�r andrew om denne. feil idenifisert?
% data(16,:) = data(10,:); % husker ikke hva dette er, hvordan putte dette
% inn i pytjon pogrammet?

% 33MeV Zinc foils:
% rows_67Cu_zn33 = [2, 3, 4, 7, 9];
% rows_62Zn_zn33 = [5, 6, 13];
% Wrongrows_62Zn_zn33 = [1, 5, 6, 12, 13];
% rows_65Ni_zn33 = [8, 17, 21];
% rows_69mZn_zn33 = [10];
% rows_63Zn_zn33 = [11, 14, 15, 20];
% rows_66Cu_zn33 = [16];
  rows_65zn_zn33 = [23];
% rows_64Cu_zn33 = [18];
% rowa_24Na_zn33 = [19, 22];

% 16MeV zirconium foils:
% rows_97Nb_zr16 = [1];
% rows_95Nb_zr16 = [5];
% rows_95Zr_zr16 = [2, 4];
% rows_97Zr_zr16 = [3];
% rows_89Zr_zr16 = [6];

% % 
% 33MeV zirconium foils:
% rows_90mY_zr33 = [1, 4]%, 8];
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
% rows_88Y_Y16 = [1, 2]%, 3];
% 
% % 33MeV Yttrium foils:
%rows_90mY_y33 = [1, 3];
%rows_87mY_y33 = [2];
%rows_87Y_y33 = [4];
%rows_88Y_y33 = [5, 6]; %7
%rows_24Na_y33 = [8];

% % 16MeV Indium foils:
% rows_116mIn_in16 = [1,3, 5, 6, 9, 10, 11, 12, 13, 14];
% rows_115mIn_in16 = [2];
% rows_114mIn_in16 = [7, 8];
% rows_113mIn_in16 = [4];
% % 
% % 33MeV Indium foils:
% rows_116mIn_In33 = [1, 5, 6, 7, 9, 11, 12, 15, 16, 17, 19, 20, 21, 22, 23, 25, 26, 27];
% rows_112mIn_In33 = [2];
% rows_111In_In33 = [3, 4];
% rows_112In_33 = [14];
% rows_115mIn_In33 = [8];
% rows_113mIn_In33 = [10];
% rows_114mIn_In33 = [13, 18];
% rows_24Na_In33 = [24, 28];

% % 16MeV Aluminum foils:
 rows_24Na_al16 = [1, 2];
% 
% % 33MeV Aluminum foils:
% rows_24Na_al33 = [1, 2];



% Select rows to plot
% varToStr = @(x) inputname(1);
rows = rows_24Na_al16;
outName = '../csv/24Na_al16MeV';
% rows = 12;
% Find rows for the desired decay product
selected_rows = data(rows,:);




% Gate on foil energy 
% 
% energy = 0; % Show all foils in one plot (not for analysis!)
% energy = 913;
% 
% % loop over all energies
% for energy = 0:0   % Show all foils in one plot (not for analysis!)
% for energy = 230:100:230   % Test mode
% for energy = 140:100:140   % Just Zirconium 16 MeV
% for energy = 240:100:240   % Just Zirconium 33 MeV
% for energy = 130:100:130   % Just Zinc 16 MeV
%for energy = 230:100:230   % Just Zinc 33 MeV
% for energy = 139:100:139   % Just Yttrium 16 MeV
%for energy = 239:100:239   % Just Yttrium 33 MeV
 for energy = 113:100:113   % Just Aluminum 16 MeV
% for energy = 213:100:213   % Just Aluminum 33 MeV
% for energy = 149:100:149   % Just Indium 16 MeV
% for energy = 349:100:349   % Just Indium 33 MeV


    if energy==0
        % Return all rows for plotting
        gammas = selected_rows;
    else
        % Oly return rows which match the desired energies
        indices = find([selected_rows.Mass] == energy);
        gammas =  selected_rows(indices);
    end
    
    % Pull the data from the GammaCounts objects
    activities     = [gammas.Activity];
    delta_ts       = [gammas.TimeSinceEoB];
    unc_activities = [gammas.UncertaintyActivity];
    
    decays     = [gammas.NumberofDecays];
    live_time  = [gammas.LiveTime];
    unc_decays = [gammas.UncertaintyNumberofDecays];
    
    errorbar(delta_ts,activities',unc_activities','.')
    
    % writing as a csv file, ' is transpose in matlab
    outfile =  [ delta_ts; activities; unc_activities; decays; unc_decays; live_time]';
    
    
    % Dump to csv for python / gnuplot
    % Turn this line on to write files out!
    % csvwrite([outName '_' num2str(energy) '.dat'],outfile);
end
