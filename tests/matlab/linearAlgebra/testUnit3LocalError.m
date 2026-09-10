function tests = testUnit3LocalError
%% SIGNATURE
% tests = testUnit3LocalError
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Check sphere-log coordinates, raw-vector Jacobians and singularity rejection.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB tests.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Validate tangent derivatives with MathCore finite differences.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeUnit3LocalError, ComputeFiniteDiffJacobian, Build3dOrthonormalBasis.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.charOldPath = path;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','..','..','matlab')));
end

function teardownOnce(testCase)
path(testCase.TestData.charOldPath);
end

function testCoordinatesAndJacobian(testCase)
dBase = [2;-3;4];
dUnit = dBase/norm(dBase);
[dAxis1,dAxis2] = Build3dOrthonormalBasis(dUnit);
for dAngle = [0,1e-7,1e-3,.4,1.7,pi-1e-3]
    dOther = 7*(cos(dAngle)*dUnit + sin(dAngle)*dAxis1);
    [dLocal,dJacobian,dBasis] = ComputeUnit3LocalError(dBase,dOther);
    verifyEqual(testCase,dLocal,[dAngle;0],'AbsTol',2e-10);
    verifyEqual(testCase,dBasis,[dAxis1,dAxis2],'AbsTol',1e-14);
    dNumeric = ComputeFiniteDiffJacobian(@(dVector) ComputeUnit3LocalError(dBase,dVector), ...
        dOther,1e-6);
    verifyEqual(testCase,dJacobian,dNumeric,'AbsTol',2e-5,'RelTol',2e-6);
    verifyEqual(testCase,dJacobian*dOther,zeros(2,1),'AbsTol',2e-10);
end
end

function testIndependentGeodesicLength(testCase)
dBase = [1;2;3];
dOther = [-4;2;1];
dLocal = ComputeUnit3LocalError(dBase,dOther);
dAngle = atan2(norm(cross(dBase,dOther)),dot(dBase,dOther));
verifyEqual(testCase,norm(dLocal),dAngle,'AbsTol',1e-14);
end

function testRejectUndefinedChart(testCase)
verifyError(testCase,@() ComputeUnit3LocalError(zeros(3,1),[1;0;0]), ...
    'ComputeUnit3LocalError:ZeroVector');
verifyError(testCase,@() ComputeUnit3LocalError([1;0;0],zeros(3,1)), ...
    'ComputeUnit3LocalError:ZeroVector');
verifyError(testCase,@() ComputeUnit3LocalError([1;0;0],[-1;0;0]), ...
    'ComputeUnit3LocalError:Antipodal');
end
