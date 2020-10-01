function traingp(xTrain, yTrain)
% Training online GP

global net gpopt ep;

% FIXED BV set - INITIALISATION
ogpemptybv(xTrain);

if gpopt.covopt.hyopt == 1
for i = 1:gpopt.postopt.nem
    ogpreset;  % remove the BV
    ogptrain(xTrain,yTrain);
end
end
ogpreset;  % remove the BV
ogppost(xTrain,yTrain);   % posterior distribution