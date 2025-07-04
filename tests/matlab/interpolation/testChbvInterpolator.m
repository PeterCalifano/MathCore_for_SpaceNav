close all
clear
clc

% TEST SETUP
charKernelsPATH = '/home/peterc/devDir/nav-backend/simulationCodes/data/SPICE_kernels';
CSPICEkerLoader(charKernelsPATH, EnumScenarioName.Itokawa);

POSVEC_POLYDEG_GT = 26;
ATTQUAT_POLYDEG_GT = 25;
dEvalPointID = 150;

projectDir = pwd;

cd(fullfile(charKernelsPATH, "common"));
cspice_furnsh('mkcommon.mk');
cd(projectDir);

% Construct data (Sun ephemerides and Body attitude)
dTimestamps = 1:10000;
dDomainLB = min(dTimestamps);
dDomainUB = max(dTimestamps);

ET0 = cspice_str2et('Apr 25 2006 06:00:00'); % From Alban
ET_SPAN = ET0 + dTimestamps;

TARGET_NAME = 'ITOKAWA';
FixedFrame = 'Itokawa_fixed'; % Check corresponding tf file
InertialFrame = 'ECLIPJ2000';

% Body attitude data
dDCM_INfromTB_seq = cspice_pxform(FixedFrame, InertialFrame, ET_SPAN);
dQuat_INfromTB_seq = DCM2quatSeq(dDCM_INfromTB_seq, false);

% Sun position
dSunPosition_IN = 1000 * cspice_spkpos('SUN', ET_SPAN, 'ECLIPJ2000', 'none', TARGET_NAME);

%% test_ChbvInterpolator
% Position vector
[dChbvCoeffs, ~, strfitStats] = fitChbvPolynomials(POSVEC_POLYDEG_GT, dTimestamps, ...
    dSunPosition_IN, dDomainLB, dDomainUB, true);

cellRefInterpolant    = cell(2,1);
cellRefInterpolant{1} = dChbvCoeffs;
cellRefInterpolant{2} = strfitStats;

% Define interpolant (MATLAB)
objSunPositionInterp = CChbvInterpolator(dTimestamps, POSVEC_POLYDEG_GT, EnumInterpType.VECTOR);

% Test polynomial terms evaluation
[objSunPositionInterp, evaluatedPoly] = objSunPositionInterp.evalPoly(dTimestamps(dEvalPointID), true);

% Test fitting
[objSunPositionInterp, dInterpCoeffsMatrix, strFitStats] = objSunPositionInterp.fitDataMatrix(dSunPosition_IN);

cellInterpolantMATLAB    = cell(2,1);
cellInterpolantMATLAB{1} = dInterpCoeffsMatrix(:);
cellInterpolantMATLAB{2} = strFitStats;

% Equivalence check
sum(abs(cellInterpolantMATLAB{1} - cellRefInterpolant{1}))

return

%% test_ChbvInterpolator_attitude

% Reference (validated functions)
% Attitude quaternion
[dChbvCoeffs, dScaledInterpDomain, dswitchIntervals, ...
    strfitStats] = fitAttQuatChbvPolynmials( ...
    ATTQUAT_POLYDEG_GT, ...
    dTimestamps, ...
    dQuat_INfromTB_seq, ...
    dDomainLB, ...
    dDomainUB, ...
    true); %#ok<*UNRCH>

cellRefAttitudeInterpolant    = cell(3,1);
cellRefAttitudeInterpolant{1} = dChbvCoeffs;
cellRefAttitudeInterpolant{2} = dswitchIntervals;
cellRefAttitudeInterpolant{3} = strfitStats;

% Define interpolant (MATLAB)
objAttitudeInterp = CChbvInterpolator(dTimestamps, ATTQUAT_POLYDEG_GT, EnumInterpType.QUAT);

% Test polynomial terms evaluation
[objAttitudeInterp, evaluatedAttPoly] = objAttitudeInterp.evalPoly(dTimestamps(dEvalPointID), true);

% Test fitting
[objAttitudeInterp, dInterpCoeffsMatrix, strFitStats] = objAttitudeInterp.fitDataMatrix(dQuat_INfromTB_seq);

cellAttInterpolantMATLAB    = cell(2,1);
cellAttInterpolantMATLAB{1} = dInterpCoeffsMatrix(:);
cellAttInterpolantMATLAB{2} = strFitStats;

% Equivalence check
sum(abs(cellAttInterpolantMATLAB{1} - cellRefAttitudeInterpolant{1}))

return




