function tests = testUnit3BasisTransport
%% SIGNATURE
% tests = testUnit3BasisTransport
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Check tangent-coordinate transport, chart changes and static generated execution.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB tests.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Validate fixed-reference tangent transport.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeUnit3BasisTransport, Build3dOrthonormalBasis, RotationVectorToDCM.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function testKnownRotationAndChartChange(testCase)
dReference = [0;0;1];
dReferenceBasis = [1,0;0,1;0,0];
dRotation = RotationVectorToDCM([0;pi/3;0]);
dDirection = dRotation*dReference;
dBasis = dRotation*dReferenceBasis;
verifyEqual(testCase,ComputeUnit3BasisTransport(dReference,dReferenceBasis, ...
    dDirection,dBasis),eye(2),'AbsTol',1e-14);
dChartChange = [cos(.7),-sin(.7);sin(.7),cos(.7)];
verifyEqual(testCase,ComputeUnit3BasisTransport(dReference,dReferenceBasis, ...
    4*dDirection,dBasis*dChartChange),dChartChange','AbsTol',1e-14);
verifyError(testCase,@() ComputeUnit3BasisTransport(dReference,dReferenceBasis, ...
    -dReference,dReferenceBasis),'ComputeUnit3LocalError:Antipodal');
end

function testClosedLoopHasNoAccumulatedRotation(testCase)
dReference = [0;0;1];
dReferenceBasis = [1,0;0,1;0,0];
dDirections = [dReference,[.5;0;sqrt(.75)],[0;.5;sqrt(.75)],dReference];
dPreviousMap = eye(2);
dBias = [.003;-.002];
dInitialBias = dBias;
for ui32Index = uint32(1:4)
    [dAxis1,dAxis2] = Build3dOrthonormalBasis(dDirections(:,ui32Index));
    dMap = ComputeUnit3BasisTransport(dReference,dReferenceBasis, ...
        dDirections(:,ui32Index),[dAxis1,dAxis2]);
    dBias = (dMap*dPreviousMap')*dBias;
    verifyEqual(testCase,dMap*dMap',eye(2),'AbsTol',1e-14);
    verifyEqual(testCase,dMap'*dBias,dInitialBias,'AbsTol',1e-16);
    dPreviousMap = dMap;
end
end

function testStaticMex(testCase)
assumeFalse(testCase,isempty(which('codegen')),'MATLAB Coder is required.');
charBuildDir = tempname;
mkdir(charBuildDir);
addpath(charBuildDir);
testCase.addTeardown(@() Cleanup_(charBuildDir));
objConfig = coder.config('mex');
objConfig.EnableVariableSizing = false;
objConfig.EnableDynamicMemoryAllocation = false;
dReference = [0;0;1];
dReferenceBasis = [1,0;0,1;0,0];
cellInputs = {dReference,dReferenceBasis,dReference,dReferenceBasis};
codegen('-config',objConfig,'ComputeUnit3BasisTransport','-args',cellInputs, ...
    '-o',fullfile(charBuildDir,'Unit3BasisTransport_mex'),'-d',fullfile(charBuildDir,'build'));
for dAngle = [0,1e-8,.6,2.8]
    dRotation = RotationVectorToDCM([dAngle;0;0]);
    cellInputs{3} = dRotation*dReference;
    cellInputs{4} = dRotation*dReferenceBasis;
    verifyEqual(testCase,Unit3BasisTransport_mex(cellInputs{:}), ...
        ComputeUnit3BasisTransport(cellInputs{:}),'AbsTol',1e-13);
end
end

function Cleanup_(charBuildDir)
clear Unit3BasisTransport_mex
rmpath(charBuildDir);
rmdir(charBuildDir,'s');
end
