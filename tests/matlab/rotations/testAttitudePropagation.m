close all
clear
clc

LoadUserPathConfig; % Does nothing for now
SetupOptions;
charProjectDir = pwd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCENARIO SETUP
%% Options and input specification
enumTrajName                = EnumTrajectoryNames.SSTO_1p4; %EnumTrajectoryNames.RTO_4t1_J13p0; % RTO_4t1_J11p0_60dt, SSTO_14_60dt, SSTO_12_60dt

% Timegrid
ui32InitialTimeID           = 1;
dFramesDeltaTime            = 120;    % [s]
dTotalTimeDuration          = 12*86400;  % [s]
ui32ImgAcquisitionStepFreq  = uint32(1);
dInitialRelTimestamp = 0; 
bUseKernelInitialTimestamp = true;
bIsTimeGridRelative = true;

% Output units
bConvertData_m2mBU_beforeRender = true;
charKernelLengthUnits           = "km";
charOutputLengthUnits           = "m";

% Frames
enumWorldFrame                  = "J2000";
enumTargetFrame                 = "APOPHIS_FIXED"; %EnumFrameName_RCS1.APOPHIS_FIXED; % "IAU_EARTH"
% enumTargetFrame                 = EnumFrameName_RCS1.APOPHIS_FIXED; % "IAU_EARTH"

%%% Trajectory kernel loader input specifications
% Target point, body and frames
ui32RenderingTargetBodyID       = uint32(0); % Keep 0 if target is the main body
varTargetID                     = int32(-19920605); % -10003001
% Target frames 
varReferenceCentre          = "APOPHIS";
varTargetBodyID             = "APOPHIS";

% Kernel name and mk path
enumTrajectKernelName       = enumTrajName; % kernelFUT_h600

% charTrajKernelFolderPath    = "/home/peterc/devDir/projects-DART/data/future/phase-C/kernels/trajectories";
charTrajKernelFolderPath = "/home/peterc/devDir/projects-DART/data/rcs-1/phase-B/kernels/trajectories";

% Timescale for timegrids from kernel
charKernelTimescale  = "ET";
charUserDefTimescale = "ET";

% Additional bodies
% cellAdditionalTargetsID          = {"MOON"};
% cellAdditionalTargetNames        = {};
% bRequireAdditionalBodiesAttitude = true(1,1);
% cellAdditionalTargetsFrames      = {"IAU_MOON"};

cellAdditionalTargetsID          = {};
cellAdditionalTargetNames        = {};
bRequireAdditionalBodiesAttitude = false(1,1);
cellAdditionalTargetsFrames      = {};

dTimegrid = 0.0:120:8*86400;
bIsTimeGridRelative = true;

charKernelPathRoot = "/home/peterc/devDir/projects-DART/data/rcs-1/phase-B/kernels/mk";

if not(strcmpi(charKernelPathRoot, ""))
    charCurrentDir = pwd;
    cd(charKernelPathRoot)
    cspice_furnsh('kernels.mk')
    cd(charCurrentDir);
end

% Load dataset data from kernels
[objDataset, ~, ~] = LoadReferenceDataFromKernels(varTargetID, ...
                                                enumTrajectKernelName, ...
                                                dTimegrid, ...
                                                enumWorldFrame, ...
                                                varReferenceCentre, ...
                                                enumTargetFrame, ...
                                                "bLoadManoeuvres", false, ...
                                                "cellAdditionalTargetsID", cellAdditionalTargetsID, ...
                                                "cellAdditionalTargetNames", cellAdditionalTargetNames, ...
                                                "charKernelLengthUnits", charKernelLengthUnits,...
                                                "charTrajKernelFolderPath", charTrajKernelFolderPath, ...
                                                "charOutputLengthUnits", charOutputLengthUnits, ...
                                                "varTargetBodyID", varTargetBodyID, ...
                                                "bIsTimeGridRelative", bIsTimeGridRelative, ...
                                                "bAdditionalBodiesRequireAttitude", bRequireAdditionalBodiesAttitude, ...
                                                "cellAdditionalTargetFrames", cellAdditionalTargetsFrames, ...
                                                "charKernelTimescale", charKernelTimescale, ...
                                                "charUserDefTimescale", charUserDefTimescale, ...
                                                "bUseKernelInitialTimestamp", bUseKernelInitialTimestamp);

% Generate new attitude data for the target, starting from initial condition
dAbsIntergTimegrid = objDataset.dTimestamps;
dQuat0      = DCM2quat(objDataset.dDCM_TBfromW(:,:,1), false);
% i32CKHandle = cspice_cklpf( char(charCKKernelName) );

% Documentation: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/mice/cspice_ckgpav.html
i32ApophisFrame_ID = -200999422; % int32 ID of the frame 
i32ApophisFrame_ID = 20099942;

[dDCM, dRefAngVelSeq, dTimeOuts, bFound] = cspice_ckgpav(int32(i32ApophisFrame_ID), ...
                                                           dAbsIntergTimegrid, ...
                                                           1E-6, ...
                                                           'APOPHIS_HF'); % Requires CK kernels!

%% Integrate using constant angular velocity
% Compute perturbed angular velocity at initial time instant
dPertubAngVel0 = dRefAngVelSeq(:,1) + deg2rad([0.01; 0.01; 0.01]) .* randn(3,1);

% Integrate attitude kinematics 
objQuatIntegr = CQuatKinematicsIntegrator();
dTimestep = 1;

[dQuatEndConstant0, dTimegridOut, dQuatSeqConstant0] = objQuatIntegr.integrate(dQuat0, ...
                                                                dPertubAngVel0, ...
                                                                dTimegrid, ...
                                                                'rkmk4', ...
                                                                dTimestep);
% Plot visualization
objDataset.plotDatasetData("dTargetAttitudeSet2", QuatSeq2DCM(dQuatSeqConstant0, false), ...
                            "bPlotTargetAttitude", true);
return

%% Integrate using angular velocity profile (internally defined spline)
% Integrate attitude kinematics 
objQuatIntegr = CQuatKinematicsIntegrator(); %#ok<UNRCH>
dTimestep = 1;

[dQuatEndConstant0, dTimegridOut, dQuatSeqConstant0] = objQuatIntegr.integrate(dQuat0, ...
                                                                dRefAngVelSeq, ...
                                                                dTimegrid, ...
                                                                'rkmk4', ...
                                                                dTimestep, ...
                                                                1.0, ...
                                                                dAbsIntergTimegrid);


% Plot visualization
objDataset.plotDatasetData("dTargetAttitudeSet2", QuatSeq2DCM(dQuatSeqConstant0, false), ...
                            "bPlotTargetAttitude", true);
return

%% Integrate dynamics and kinematics
% TODO


