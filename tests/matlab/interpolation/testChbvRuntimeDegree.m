function tests = testChbvRuntimeDegree
%% SIGNATURE
% tests = testChbvRuntimeDegree
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Verify active Chebyshev degrees within fixed storage, coefficient packing and static MEX execution.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB tests.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Cover runtime degrees and fixed maximum allocation.
% 10-09-2026  Pietro Califano, Codex gpt-6    Check quaternion series at changing runtime degrees.
% 11-09-2026  Pietro Califano, Codex gpt-6    Remove unused runtime sign-switch metadata.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% evalChbvPolyWithCoeffs, EvalRecursiveChbv, MATLAB Coder for the MEX test.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function testActiveDegreeBelowCapacity(testCase)
for ui32Degree = uint32([2,3,8])
    [dStorage,dActive] = MakeCoefficients_(ui32Degree);
    for dTime = [-2,-.3,1,4]
        dScaledTime = (2*dTime-2)/6;
        dBasis = cos((0:double(ui32Degree))'*acos(dScaledTime));
        dExpected = reshape(dActive,double(ui32Degree)+1,3)'*dBasis;
        dActual = evalChbvPolyWithCoeffs(ui32Degree,uint32(3),dTime,dStorage,-2,4, ...
            uint32(numel(dActive)),uint32(8));
        verifyEqual(testCase,dActual,dExpected,'AbsTol',2e-13);
        verifyEqual(testCase,evalChbvPolyWithCoeffs(ui32Degree,uint32(3), ...
            dTime,dActive,-2,4),dExpected,'AbsTol',2e-13);
    end
end
end

function testRejectInvalidCapacity(testCase)
verifyError(testCase,@() evalChbvPolyWithCoeffs(uint32(9),uint32(3),0,zeros(30,1), ...
    -2,4,uint32(30),uint32(8)), 'evalChbvPolyWithCoeffs:DegreeExceedsMaximum');
verifyError(testCase,@() evalChbvPolyWithCoeffs(uint32(3),uint32(3),0,zeros(9,1), ...
    -2,4,uint32(12),uint32(8)), 'evalChbvPolyWithCoeffs:CoefficientCount');
verifyError(testCase,@() EvalRecursiveChbv(uint32(9),0,uint32(8)), ...
    'EvalRecursiveChbv:DegreeExceedsMaximum');
end

function testQuaternionRuntimeDegree(testCase)
% Compare against the explicit cosine basis; unused capacity must not affect the output.
for ui32Degree = uint32(2:8)
    dActive = sin((1:4 * (double(ui32Degree) + 1))');
    dStorage = nan(36, 1);
    dStorage(1:numel(dActive)) = dActive;
    for dTime = [-2, -0.3, 4]
        dBasis = cos((0:double(ui32Degree))' * acos((2 * dTime - 2) / 6));
        dExpected = reshape(dActive, double(ui32Degree) + 1, 4)' * dBasis;
        dActual = evalAttQuatChbvPolyWithCoeffs(ui32Degree, uint32(4), dTime, ...
            dStorage, -2, 4, uint32(8));
        verifyEqual(testCase, dActual, dExpected, 'AbsTol', 2e-13);
    end
end
end

function testQuaternionStaticMex(testCase)
assumeFalse(testCase, isempty(which('codegen')), 'MATLAB Coder is required.');
charBuildDir = tempname;
mkdir(charBuildDir);
addpath(charBuildDir);
testCase.addTeardown(@() CleanupMex_(charBuildDir));
objConfig = coder.config('mex');
objConfig.EnableVariableSizing = false;
objConfig.EnableDynamicMemoryAllocation = false;
cellInputs = {uint32(2), coder.Constant(uint32(4)), 0, zeros(36, 1), ...
    -2, 4, coder.Constant(uint32(8))};
codegen('-config', objConfig, 'evalAttQuatChbvPolyWithCoeffs', '-args', cellInputs, ...
    '-o', fullfile(charBuildDir, 'QuatRuntimeDegree_mex'), '-d', fullfile(charBuildDir, 'build'));

% Reuse one binary with different active degrees and NaNs in the unused tail.
for ui32Degree = uint32(2:8)
    dActive = sin((1:4 * (double(ui32Degree) + 1))');
    dStorage = nan(36, 1);
    dStorage(1:numel(dActive)) = dActive;
    for dTime = [-2, -0.3, 4]
        dBasis = cos((0:double(ui32Degree))' * acos((2 * dTime - 2) / 6));
        dExpected = reshape(dActive, double(ui32Degree) + 1, 4)' * dBasis;
        dActual = QuatRuntimeDegree_mex(ui32Degree, uint32(4), dTime, ...
            dStorage, -2, 4, uint32(8));
        verifyEqual(testCase, dActual, dExpected, 'AbsTol', 2e-13);
    end
end
verifyError(testCase, @() QuatRuntimeDegree_mex(uint32(9), uint32(4), 0, ...
    zeros(36, 1), -2, 4, uint32(8)), ...
    'evalChbvPolyWithCoeffs:DegreeExceedsMaximum');
end

function testStaticMexWithRuntimeDegree(testCase)
assumeFalse(testCase,isempty(which('codegen')),'MATLAB Coder is required.');
charBuildDir = tempname;
mkdir(charBuildDir);
addpath(charBuildDir);
testCase.addTeardown(@() CleanupMex_(charBuildDir));
objConfig = coder.config('mex');
objConfig.EnableVariableSizing = false;
objConfig.EnableDynamicMemoryAllocation = false;
cellInputs = {uint32(2),coder.Constant(uint32(3)),0,zeros(27,1),-2,4, ...
    uint32(9),coder.Constant(uint32(8))};
codegen('-config',objConfig,'evalChbvPolyWithCoeffs','-args',cellInputs, ...
    '-o',fullfile(charBuildDir,'ChbvRuntimeDegree_mex'),'-d',fullfile(charBuildDir,'build'));
for ui32Degree = uint32(2:8)
    [dStorage,dActive] = MakeCoefficients_(ui32Degree);
    for dTime = [-2,-.3,4]
        dExpected = evalChbvPolyWithCoeffs(ui32Degree,uint32(3),dTime,dActive,-2,4);
        dActual = ChbvRuntimeDegree_mex(ui32Degree,uint32(3),dTime,dStorage,-2,4, ...
            uint32(numel(dActive)),uint32(8));
        verifyEqual(testCase,dActual,dExpected,'AbsTol',2e-13);
    end
end
verifyError(testCase,@() ChbvRuntimeDegree_mex(uint32(9),uint32(3),0,zeros(27,1), ...
    -2,4,uint32(30),uint32(8)), 'evalChbvPolyWithCoeffs:DegreeExceedsMaximum');
end

function [dStorage,dActive] = MakeCoefficients_(ui32Degree)
dActive = sin((1:3*(double(ui32Degree)+1))');
% Only the packed active prefix is meaningful; unused capacity must never be read.
dStorage = nan(27,1);
dStorage(1:numel(dActive)) = dActive;
end

function CleanupMex_(charBuildDir)
clear ChbvRuntimeDegree_mex
clear QuatRuntimeDegree_mex
rmpath(charBuildDir);
rmdir(charBuildDir,'s');
end
