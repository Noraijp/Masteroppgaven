    activities     = [obj.Activity];
    delta_ts       = [obj.TimeSinceEoB];
    unc_activities = [obj.UncertaintyActivity];
    
%     errorbar(delta_ts,activities',unc_activities','.')
    
    % writing as a csv file, ' is transpose in matlab
    outfile =  [ delta_ts; activities; unc_activities]';
    
    
    % Dump to csv for python / gnuplot
    % Turn this line on to write files out!
    csvwrite(['../csv/65Zn_zn16MeV_' num2str(energy) '.dat'],outfile);