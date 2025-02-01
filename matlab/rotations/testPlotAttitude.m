close all
clear
clc


%% Plotting routines test script

addpath(genpath("S:/TestInputData"));
addpath(genpath("./TestScripts"));
addpath(genpath("./CREmodules"));
addpath(genpath("./AnalyticalUKF"))

scenarioData = 'exp';
[dataPALT, dataDKE, dataTNAV, dataADCS, dataGNC, dataIP, dataFDIR, dataEPH] = loadTestData(scenarioData);
plotData = 0; % Enables plot of state wrt D2 (RefData)

% Load data
getParamsRefData;

%% TEST
testNumber = 0;

switch testNumber
    case 0

        Nlast = 120000;
        % Plot Attitude Quaternion
%         [fig, o_dDCM_Target2Fixed] = plotAttitudeQuat(i_dQuatSeq, i_dOriginPos, i_bConvFlag, i_bPlotFrame)
        i_dOriginPos = R_SC_True(:, 1:Nlast);
%         i_dOriginPos = zeros(3, 1);
        i_bConvFlag = 0;
        i_bPlotFrame = false;
        i_dQuatSeq = qCAMwrtIN_ref(:, 1:Nlast);

        plotAttitudeQuat(i_dQuatSeq, i_dOriginPos, i_bConvFlag, i_bPlotFrame);

    case 1



end