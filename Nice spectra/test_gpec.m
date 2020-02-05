clear all; clc;

filename = 'AI081618_10cm_Zn02.Spe';

file = importdata(filename,'',12);
spectrum = file.data;

% Give your coefficients for energy calibration here
% 
% Offset is 0th-order (constant) term
% Gain is 1st-order (linear) term
% 
offset = -2.7512E-02;
gain   = 1.8736E-001;
En = offset +   (gain.*(0:length(spectrum)-1));

h1 = semilogy(En,spectrum, 'Color', [179, 0, 0]./255);
set(get(get(h1,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend

% Manually set axis limits, if desired
% 
xlim([0 1500])    % in keV
% ylim([1e0 1e5])

% Find "strong" peaks
[pks,locs] =  findpeaks(spectrum,'MinPeakProminence',50);

% pks = pks([21,24,25]);
% locs = locs([21,24,25]);

% Energy calibration for peak locations
locs = offset +   (gain.*locs);
% i = locs<377;
% locs(i) = [];
% pks(i) = []; 
%  Remove 511 keV peak
i = (locs>508 & locs<513);
% i = [1,0,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0]==1;
locs(i) = [];
pks(i) = []; 



% 
% % % Manual addition of additional peaks
% % locs(30) = 263;
% % pks(30) = 1692;
% % locs(31) = 1471;
% % pks(31) = 102;

locs(22) = 596.7;
pks(22) = 23;

locs(23) = 1345;
pks(23) = 22;
% 
% locs(8) = 1314;
% pks(8) = 15;
% 
% locs(9) = 1366;
% pks(9) = 25;
% 
% locs(10) = 1406;
% pks(10) = 10;
% 
% 
hold on

 h2 = semilogy(locs([5, 22]),pks([5, 22]).*2,'s');
    set(h2, 'MarkerFaceColor', get(h2, 'Color'));
%     '$^{87m}$Y IT Decay',...
             
      h4 = semilogy(locs([11 12 16 18]),pks([11 12 16 18]).*2,'d');
    set(h4, 'MarkerFaceColor', get(h4, 'Color'));
%     '$^{88}$Y $\epsilon$ Decay',...
             
    h6 = semilogy(locs(19),pks(19).*2.2,'p');
    set(h6, 'MarkerFaceColor', get(h6, 'Color'));
%     '$^{89}$Zr $\epsilon$ Decay',...
             
     h7 = semilogy(locs([20]),pks([20]).*1.7,'o');
    set(h7, 'MarkerFaceColor', get(h7, 'Color'));
%     '$^{90}$Mo $\epsilon$ Decay',...
             
% %      h8 = semilogy(locs([2,3,8,14,15,19,21,23,24,25,26,27]),pks([2,3,8,14,15,19,21,23,24,25,26,27]).*2,'*', 'linewidth',1.2);
% %      h8 = semilogy(locs([2,3,8,15,19,21,22,23,24,26,27]),pks([2,3,8,15,19,21,22,23,24,26,27]).*2.3,'*', 'linewidth',1.2);
%    
    h8 = semilogy(locs(23),pks(23).*2.3,'*', 'linewidth',1.2);
    set(h8, 'MarkerFaceColor', get(h8, 'Color'));
% '$^{90}$Nb $\epsilon$ Decay',...
             
    h9 = semilogy(locs(21),pks(21).*1.3,'+', 'linewidth',2);
    set(h9, 'MarkerFaceColor', get(h9, 'Color'));
% %     '$^{92m}$Nb $\epsilon$ Decay',...
%              
% % %      h11 = semilogy(locs([13,23]),pks([13,23]).*2,'h');
% %      h11 = semilogy(locs([13,30,31]),pks([13,30,31]).*2,'^');
% %     set(h11, 'MarkerFaceColor', get(h11, 'Color'));
% %     '$^{93m}$Mo $\epsilon$ Decay',..


%      h7 = semilogy(locs,pks.*1.7,'o');
%     set(h7, 'MarkerFaceColor', get(h7, 'Color'));
% %     Plotting all peaks

pbaspect([2 1 1])
hy = ylabel('Counts', 'FontSize', 20);
hx = xlabel('Energy (keV)', 'FontSize', 20);
hold off


hl1 = legend('$^{62}$Zn $\epsilon$ Decay',...
             '$^{67}$Cu $\beta^-$ Decay',...%              '$^{88}$Y $\epsilon$ Decay',...
             '$^{69m}$Zn $IT$ Decay',...
             '$^{65}$Zn  $\epsilon$ Decay',...
             '$^{64}$Cu $\epsilon$ Decay',...
             '$^{24}$Na $\beta^-$ Decay',...
             'location','northeast');

set([hx,hy,hl1],'Interpreter','latex')

% set(gcf,'units','points','position',[1239 116 560 420])
set(gca,'XMinorTick','on')

% 
% export_fig ./33MeV_longcount_Zn.pdf -pdf -transparent
% 
% print('./33MeV_longcount_Zn', '-dpng', '-r300');





%%


clear all; clc;

filename = 'AB0811118_10cm_Zn01.txt';

file = importdata(filename);
spectrum = file(:,2);

% Give your coefficients for energy calibration here
% 
% Offset is 0th-order (constant) term
% Gain is 1st-order (linear) term
% 
offset = -2.7512E-02;
gain   = 1.8736E-001;
En = offset +   (gain.*(0:length(spectrum)-1));

h1 = semilogy(En,spectrum, 'Color', [179, 0, 0]./255);
set(get(get(h1,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend

% Manually set axis limits, if desired
% 
xlim([0 1500])    % in keV
% ylim([1e0 1e5])

% Find "strong" peaks
[pks,locs] =  findpeaks(spectrum,'MinPeakProminence',10);

% pks = pks([21,24,25]);
% locs = locs([21,24,25]);

% Energy calibration for peak locations
locs = offset +   (gain.*locs);
% i = locs<377;
% locs(i) = [];
% pks(i) = []; 
%  Remove 511 keV peak
i = (locs>508 & locs<513);
% i = [1,0,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0]==1;
locs(i) = [];
pks(i) = []; 



% 
% % % Manual addition of additional peaks
% % locs(30) = 263;
% % pks(30) = 1692;
% % locs(31) = 1471;
% % pks(31) = 102;
% 
% locs(22) = 596.7;
% pks(22) = 23;
% 
% locs(23) = 1345;
% pks(23) = 22;
% 
% locs(8) = 1314;
% pks(8) = 15;
% 
% locs(9) = 1366;
% pks(9) = 25;
% 
% locs(10) = 1406;
% pks(10) = 10;
% 
% 
hold on

 h2 = semilogy(locs([101 103 224]),pks([101 103 224]).*2,'s');
    set(h2, 'MarkerFaceColor', get(h2, 'Color'));
% %     '$^{87m}$Y IT Decay',...
             
      h4 = semilogy(locs(544),pks(544).*2,'d');
    set(h4, 'MarkerFaceColor', get(h4, 'Color'));
% %     '$^{88}$Y $\epsilon$ Decay',...
%              
% %    
             
     h7 = semilogy(locs(603),pks(603).*1.7,'o');
    set(h7, 'MarkerFaceColor', get(h7, 'Color'));
% % %     '$^{90}$Mo $\epsilon$ Decay',...
    
    h6 = semilogy(locs(613),pks(613).*2.2,'p');
    set(h6, 'MarkerFaceColor', get(h6, 'Color'));
% % %     '$^{89}$Zr $\epsilon$ Decay',...
             
% %      h8 = semilogy(locs([2,3,8,14,15,19,21,23,24,25,26,27]),pks([2,3,8,14,15,19,21,23,24,25,26,27]).*2,'*', 'linewidth',1.2);
% %      h8 = semilogy(locs([2,3,8,15,19,21,22,23,24,26,27]),pks([2,3,8,15,19,21,22,23,24,26,27]).*2.3,'*', 'linewidth',1.2);
%    
    h8 = semilogy(locs(616),pks(616).*2.3,'*', 'linewidth',1.2);
    set(h8, 'MarkerFaceColor', get(h8, 'Color'));
% % '$^{90}$Nb $\epsilon$ Decay',...
             
    h9 = semilogy(locs(619),pks(619).*1.3,'+', 'linewidth',2);
    set(h9, 'MarkerFaceColor', get(h9, 'Color'));
% % %     '$^{92m}$Nb $\epsilon$ Decay',...
%              
% % %      h11 = semilogy(locs([13,23]),pks([13,23]).*2,'h');
% %      h11 = semilogy(locs([13,30,31]),pks([13,30,31]).*2,'^');
% %     set(h11, 'MarkerFaceColor', get(h11, 'Color'));
% %     '$^{93m}$Mo $\epsilon$ Decay',..

% 
%      h7 = semilogy(locs,pks.*1.7,'o');
%     set(h7, 'MarkerFaceColor', get(h7, 'Color'));
%     Plotting all peaks

pbaspect([2 1 1])
hy = ylabel('Counts', 'FontSize', 20);
hx = xlabel('Energy (keV)', 'FontSize', 20);
hold off


hl1 = legend('$^{67}$Cu $\beta^-$ Decay',...
             '$^{69m}$Zn $IT$ Decay',...%              '$^{88}$Y $\epsilon$ Decay',...
             '$^{63}$Zn $\epsilon$ Decay',...
             '$^{65}$Ni  $\beta^-$ Decay',...
             '$^{64}$Cu $\epsilon$ Decay',...
             '$^{24}$Na $\beta^-$ Decay',...
             'location','northeast');

set([hx,hy,hl1],'Interpreter','latex')

% set(gcf,'units','points','position',[1239 116 560 420])
set(gca,'XMinorTick','on')


export_fig ./16MeV_longcount_Zn.pdf -pdf -transparent

print('./16MeV_longcount_Zn', '-dpng', '-r300');








