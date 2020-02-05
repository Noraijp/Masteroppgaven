clear all; clc;
fclose('all');

% load effi_mat.mat

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

% Old al_key_energies
% al_key_energies = 1274.51;

% Load file names

%zn_mat_16MeV = zeros(15,25);


%
% Change element prefix for energy identification, and pre-allocate size
%
%
zn_fn = textscan(fopen('zn_fnames.txt'),'%s');   
zr_fn = textscan(fopen('zr_fnames.txt'),'%s');   
y_fn = textscan(fopen('y_fnames.txt'),'%s');   
al_fn = textscan(fopen('al_fnames.txt'),'%s');   
in_fn = textscan(fopen('in_fnames.txt'),'%s');   



% Pre-allocate size of arrays
% 33 MeV
zn_energies = zeros(2000,1);
zr_energies = zeros(2000,1);
y_energies = zeros(2000,1);
al_energies = zeros(2000,1);
in_energies = zeros(2000,1);

%             1XX        2XX      3XX
zr_rhodrs = [0.636     0.635];
zn_rhodrs = [0.663     0.646];
al_rhodrs = [0.197     0.196];
in_rhodrs = [0.358     0.000     0.354];
y_rhodrs =  [0.402     0.426];



% select which foils to analyze
fn = in_fn;
energies = in_energies;
% rhodrs = zn_rhodrs;
% mu_attenuation = load('../zn_xcom.txt');



%Fit attenuation data to smooth function

% pp = spline(mu_attenuation(:,1).*1e3,mu_attenuation(:,2));
% pp2 = pchip(mu_attenuation(:,1).*1e3,mu_attenuation(:,2));
% xp=logspace(0,4,10000);
% yp=ppval(pp,xp);
% yp2=ppval(pp2,xp);
% loglog(mu_attenuation(:,1).*1e3,mu_attenuation(:,2),'--',xp,yp2);
% xlabel('Gamma-Ray Energy (keV)')
% ylabel('\mu (cm^2/g)')
%clf;



%
% End of selection
%



% Plotting extracted unique peak energies



% Changing indices for vector insertion
a = 1;
b = 1;

for i=1:length(fn{1})
% for i=1:1
    % for i=1:3
    fname = fn{1}{i,1};
    fid = fopen(fname);
    
%     C2 = textscan(fid,'%d %f %f %f %f %d %f %f %f %f %f %s','headerlines',30);
    C2 = textscan(fid,'%f %f %f %f %f %d %f %f %f %d','headerlines',30);
    b = a + length(C2{1,1}) - 1;
    
    energies(a:b) =  C2{1,1};
    a = a + length(C2{1,1});
    
    fclose(fid);
end

energies = sort(energies);

%
% Be careful - veto  peaks
%
% energies(913) = [];
% energies(1438) = [];
% energies(773) = [];
%
%  End being careful
%

xp = 1:length(energies);
% plot(1:length(nb_energies),nb_energies,'.')
%
% Set tolerance (in keV) for peak uniqueness
%
tol = 1.5;

[C,IA] = uniquetol(energies, tol , 'DataScale', 1,'OutputAllIndices',true);


hold on
for k = 1:length(IA)
    plot(xp(IA{k}), energies(IA{k}), '.')
    meanAi = mean(energies(IA{k},:));
    C(k) = meanAi;
    plot(find((energies-meanAi)>=0,1), meanAi, 'xr')
end

xlabel('Peak Index')
ylabel('Gamma-Ray Energy (keV)')
unique_energies = C;

% close all;


% nb_fn = string(nb_fn{1})
% Extract name i and remove .dat extension
% nb_fn{1}{i,1}(1:length(nb_fn{1}{i,1})-4)


% % Convert to array of strings
% <change  char()  to string()  after you update matlab!!!
zr_fn = char(zr_fn{1});
in_fn = char(in_fn{1});
al_fn = char(al_fn{1});
zn_fn = char(zn_fn{1});
y_fn = char(y_fn{1});




%%

% unique shelves:  1     3     5    10    14    15    18    22    30    40    60
% 
%             B0               B1              B2                B3                 B4
B_5cm  = [4.09797847e+01   4.23032792e+00   1.02954775e-01   8.35777814e-03   1.19335648e+00];
B_10cm = [5.86527950e+06   1.62999681e+01   3.93652148e-02   5.60461201e-03   1.22575525e+00];
B_18cm = [1.07096818e+28   6.58778946e+01   1.14064685e-02   3.48029377e-03   1.30117611e+00];
B_22cm = [2.68177003e+72   1.68299183e+02   4.69514974e-03   4.08167495e-03   1.25650245e+00];

eff_5  = @(x)(B_5cm(1).*  exp(-B_5cm(2).*(x.^B_5cm(3)))   .*  (1-exp(-B_5cm(4).*(x.^B_5cm(5))))  )  ;
eff_10 = @(x)(B_10cm(1).*  exp(-B_10cm(2).*(x.^B_10cm(3)))   .*  (1-exp(-B_10cm(4).*(x.^B_10cm(5))))  )  ;
eff_18 = @(x)(B_18cm(1).*  exp(-B_18cm(2).*(x.^B_18cm(3)))   .*  (1-exp(-B_18cm(4).*(x.^B_18cm(5))))  )  ;
eff_22 = @(x)(B_22cm(1).*  exp(-B_22cm(2).*(x.^B_22cm(3)))   .*  (1-exp(-B_22cm(4).*(x.^B_22cm(5))))  )  ;

% beta=[ 0.1613   -2.8769    7.6767];
% eff_3n = @(x)(exp(beta(1).*(log(x)).^2 + beta(2).*log(x) + beta(3)))  ;

xp = 62:1:2500;

% semilogy(xp, eff_1(xp),xp, eff_5(xp), xp, sqrt(eff_1(xp).*eff_5(xp)), xp, eff_3n(xp))
% semilogy( xp, sqrt(eff_1(xp).*eff_5(xp)), xp, eff_3n(xp))

% plot(xp, sqrt(eff_1(xp).*eff_5(xp)) -  eff_3n(xp))
% plot(xp, (sqrt(eff_1(xp).*eff_5(xp)) -  eff_3n(xp))./eff_3n(xp))

semilogy(xp, eff_5(xp),xp, eff_10(xp),xp, eff_18(xp), xp, eff_22(xp))
legend('5cm', '10cm', '18cm','22cm')

% ylim([6E-5 3E-1])

%%

clc;


% Select foils to analyze
fn = in_fn;
key_energies = in_key_energies;
glines = in_glines;
%glines = zn_key_energies1; for in_glines
rhodrs = in_rhodrs;
mu_attenuation = load('../in_xcom.txt');

pp2 = pchip(mu_attenuation(:,1).*1e3,mu_attenuation(:,2)); %peacewise polynomial

%1e3 = grams to miligrams 


number_of_files = size(fn);
number_of_files = number_of_files(1);

% pre-allocate size of data matrix
data_mat = zeros(length(key_energies), 1+ (8 * number_of_files));
data_mat(:,1) = key_energies;



% shelves = zeros(1,length(fn));
shelves = zeros(1,number_of_files);
% detectors = shelves;



% for i=1:length(fn)
for i=1:number_of_files
    
%     fname = char(fn(i,:)); % fn = filenames in all filenemes
    fname = fn(i,:); 
    fid = fopen(strtrim(fname));
    
    % Get raw text for header regex extraction
    raw_str = fileread(strtrim(fname)); 
    
    % Parse column data to cell structure
    %   '37' refers to the number of header lines in the report file 
    C3 = textscan(fid,'%f %f %f %f %f %d %f %f %f %d','headerlines',30); %textscan, the giant block of the data for the current file.
    %C3 is a cellarray, C3{1,1} = 
    fclose(fid);
    
    % Get shelf position, for efficiency correction
    shelf = regexp(raw_str, 'Shelf:\s+(\d+)[on]?', 'tokens');   % finding shelf number form the report files, regex= regular exression(a type of search)
    shelves(i) = cell2mat(cellfun(@(x) str2double(x{:}), shelf, 'UniformOutput', false)); % fTakes the rest of the text og that line after the coloumn and converts in from a string to a number
    shelf = shelves(i);  % Taking the number and saving it to a veriable called shelfs
    
%     % Get detector ID, for efficiency correction
%     detector = regexp(raw_str, 'Detector:\s+(\d+)[on]?', 'tokens');
%     detectors(i) = cell2mat(cellfun(@(x) str2double(x{:}), detector, 'UniformOutput', false));
%     detector = detectors(i);
    
%     Select efficiency curve based on shelf
%     if detector==1
        if shelf==5
            %effcal = eff_5;
            effcal = 'eff_room131_5.mat';
        elseif shelf==10
            %effcal = eff_10;
            effcal = 'eff_room131_10.mat';
        elseif shelf==18
            %effcal = eff_18;
            effcal = 'eff_room131_18.mat';
        elseif shelf==22
            %effcal = eff_22;
            effcal = 'eff_room131_22.mat';
%         elseif shelf==14
%             effcal = eff_14;
%         elseif shelf==15
%             effcal = eff_15;
%         elseif shelf==18
%             effcal = eff_18;
%         elseif shelf==22
%             effcal = eff_22;
%         elseif shelf==30
%             effcal = eff_30;
%         elseif shelf==40
%             effcal = eff_40;
%         elseif shelf==60
%             effcal = eff_60;
        end
%     elseif detector==2
%         if shelf==1
%             effcal = eff_1;
%         end     
%     end
    
    % Turn me off!!!  Just for testing Hannah output
%     effcal=1;
    
    % Number of columns - testing shelf regex
%     num_cols = 10;
    num_cols = 9;
    
    
    extract_energies = C3{1,1};
    dim = length(C3{1,1});
    
    % Make matrix to hold data for this spectrum
    sto_mat = zeros(dim,num_cols);
    
    
    % Platform: 1 Type:
    % pattern = 'Platform:\s+(\d+)';
    % data = regexp(indata, pattern, 'tokens');
    % data = cell2mat(cellfun(@(x) str2double(x{:}), data, 'UniformOutput', false));
    
    % Extract header
    mass = regexp(raw_str, 'Mass:\s+(\d+)', 'tokens'); % findinf foil ID
    mass = cell2mat(cellfun(@(x) str2double(x{:}), mass, 'UniformOutput', false)); 
    mass_str = num2str(mass);
    foil_id = str2num(mass_str(1:end-2));
    
    ct = regexp(raw_str, 'Datetime:\s+([-\w: ]+)', 'tokens'); % find count time
%     ct{1,1}{1,1}
%     ct = cell2mat(cellfun(@(x) str2double(x{:}), ct, 'UniformOutput', false))
    ct = juliandate(datetime(ct{1,1}{1,1}));
    
    cl = regexp(raw_str, 'Live:\s+(\d+)', 'tokens'); % find livetime
    cl = cell2mat(cellfun(@(x) str2double(x{:}), cl, 'UniformOutput', false));
    
    % Append data from header to columns, putting the data from the big
    % data and putting it into a big matrix
    sto_mat(:,4) = mass;
    sto_mat(:,5) = ct;
    sto_mat(:,6) = cl;
    
    % t12 = regexp(raw_str, 'Mass:\s+(\d+)', 'tokens');
    % t12 = cell2mat(cellfun(@(x) str2double(x{:}), t12, 'UniformOutput', false))
    %
    % Igamma = regexp(raw_str, 'Mass:\s+(\d+)', 'tokens');
    % Igamma = cell2mat(cellfun(@(x) str2double(x{:}), Igamma, 'UniformOutput', false))

    
    % Calculate the correct efficiency
    [eff, unc_eff, pcov, popt] = efficiency_calibration(C3{1,1}, effcal); % c3{1,1} = E_gamma, alle gammaenergiene i rapporten

    
    % Pull out gamma-ray energies and net counts
%                     energy        net counts (efficiency and attenuation corrected)                                      %error
    sto_mat(:,1:3) = [C3{1,1}      double(C3{1,6})./(effcal(C3{1,1}).*exp(-ppval(pp2,C3{1,1}).* 0.5 .* rhodrs(foil_id)) )       C3{1,7}];
% C3{1,1} = E_gamma. 
% C3{1,7}= unc in number of counts, 
% C3{1,6} = net number of cunts (corrected for efficiency and gammarays
% self intanuttuon)
%C3{1,6} = net counts corrected for efficiency and self attinuation, then
%divide by effcal(C3{1,1}) and e^{-mu rho dr/2}
%0.5 gammarays are born half way into the foil, we assume.
%                       
%     sto_mat(:,1:3) = [C3{1,1}      double(C3{1,6})./(effcal.*exp(-ppval(pp2,C3{1,1}).* 0.5e-3 .* rhodrs(foil_id)) )       C3{1,7}];
%     sto_mat(:,10) = shelf;
    
    
    
    % Replace energy values with closest from key energies, otherwise delete?
    % Indices of rows to keep
    ind = zeros(2,dim);
    
    for j=1:dim
        [Ydiff, idx] = min(abs(key_energies - sto_mat(j,1)));
        if Ydiff <= 1.0
            ind(:,j) = [j;idx];
        end
    end
    
    % Remove null values
    % ind = ind(ind ~= 0);
    ind = ind(:,ind(1,:) ~= 0);
    
    % Extract rows within 1 keV or key energies
    sto_mat = sto_mat(ind(1,:), :);
    % Replace with key energies
    sto_mat(:,1) = key_energies(ind(2,:));
    % append the t12 and Igamma data from a key
    sto_mat(:,7:9) = glines(ind(2,:),:);
    
    
    
    % Merge matrix into data_mat
    
    % [~,ii] = ismember(A(:,2),B(:,1));
    % [~,ii] = ismember(sto_mat(:,1),key_energies(:,1));
    % out = [key_energies, sto_mat(ii,2)];
    
    % data_mat(ind(2,:),(2:7)+7.*(i-1)+1) = sto_mat;
    
    % Need to append the t12 and Igamma data from a key ?
    % i=1;
%     data_mat(ind(2,:),(1:8)+8.*(i-1)+1) = sto_mat(:,2:length(sto_mat));
    data_mat(ind(2,:),(1:8)+8.*(i-1)+1) = sto_mat(:,2:9);

    % data_mat(ind(2,:),((1:5)+5.*(i-1)+2):((1:5)+5.*(i-1)+4)) = sto_mat(:,2:length(sto_mat))
    
end

unique(shelves)



fclose('all');

