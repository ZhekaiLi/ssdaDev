disp('  Begin to run case'), disp(analysisopt.analysistype), disp(' for '), disp(analysisopt.runTimes), disp(' times');
analysisopt.echo_flag = 2;
for i = 1:analysisopt.runTimes 
    ferum_ZK;
    result.accGPRate(i) = sum(ssda_result.ssda_Data.AccRate(2,:));
    result.accSSRate(i) = sum(ssda_result.ssda_Data.AccRate(3,:));
    result.GPdividedBySS(i) = result.accGPRate(i)/result.accSSRate(i);
end