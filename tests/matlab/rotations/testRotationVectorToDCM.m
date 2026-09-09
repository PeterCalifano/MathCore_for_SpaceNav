function tests = testRotationVectorToDCM
%% SIGNATURE
% tests = testRotationVectorToDCM
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Validate the rotation exponential and its left Jacobian with independent
% matrix exponentials, integral identities and MathCore central differences.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% None.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tests    Function-based MATLAB test suite.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Validate optional left-Jacobian output.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% RotationVectorToDCM, ComputeFiniteDiffJacobian, skewSymm, ValidateDCM, DCM2quat.
% -------------------------------------------------------------------------------------------------------------
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
charTestFolder = fileparts(mfilename('fullpath'));
charMatlabFolder = fullfile(charTestFolder, '..', '..', '..', 'matlab');
charRotationsFolder = fullfile(charMatlabFolder, 'rotations');
testCase.TestData.charOriginalPath = path;
addpath(charMatlabFolder, charRotationsFolder, fullfile(charRotationsFolder, 'dcmLib'), ...
    fullfile(charRotationsFolder, 'quatLib'), fullfile(charMatlabFolder, 'linearAlgebra'), ...
    fullfile(charMatlabFolder, 'misc'));
end

function teardownOnce(testCase)
path(testCase.TestData.charOriginalPath);
end

function testZeroRotationVectorReturnsIdentity(testCase)
[dDCM, dLeftJacobian] = RotationVectorToDCM(zeros(3, 1));

verifyEqual(testCase, dDCM, eye(3), "AbsTol", 0.0);
verifyEqual(testCase, dLeftJacobian, eye(3), "AbsTol", 0.0);
verifyTrue(testCase, ValidateDCM(dDCM));
end

function testSmallAngleRotationVectorReturnsValidDCM(testCase)
dRotationVector = [1.0e-4; 2.0e-4; 3.0e-4];

dDCM = RotationVectorToDCM(dRotationVector);
dQuat = DCM2quat(dDCM, false);

verifyTrue(testCase, ValidateDCM(dDCM));
verifyLessThanOrEqual(testCase, max(abs(vecnorm(dDCM).^2 - 1.0)), eps('single'));
verifyLessThanOrEqual(testCase, norm(cross(dDCM(:, 1), dDCM(:, 2)) - dDCM(:, 3)), eps('single'));
verifyEqual(testCase, size(dQuat), [4, 1]);
verifyLessThanOrEqual(testCase, abs(norm(dQuat) - 1.0), double(eps('single')));
end

function testLeftJacobianAgainstSharedFiniteDifferences(testCase)
dAxis = [1; -2; 3]/sqrt(14);
for dAngle = [0, 1e-9, 0.999e-4, 1.001e-4, 0.2, 1.7, pi-1e-4, 4.0]
    dVector = dAngle*dAxis;
    [dDCM, dJacobian] = RotationVectorToDCM(dVector);
    dMatrixDerivative = ComputeFiniteDiffJacobian(@RotationVectorToDCM, dVector, 1e-6);
    dNumericJacobian = zeros(3);

    % Convert dR/dphi to a left local rotation via (dR/dphi)*R'.
    for ui32Column = uint32(1):uint32(3)
        dLocalSkew = reshape(dMatrixDerivative(:,ui32Column),3,3)*dDCM';
        dNumericJacobian(:,ui32Column) = ...
            [dLocalSkew(3,2); dLocalSkew(1,3); dLocalSkew(2,1)];
    end

    verifyEqual(testCase, dJacobian, dNumericJacobian, 'AbsTol', 2e-9);
    verifyEqual(testCase, dDCM, expm(skewSymm(dVector)), 'AbsTol', 2e-14);
    verifyEqual(testCase, RotationVectorToDCM(dVector), dDCM);
end
end

function testLeftJacobianMatchesIntegralIdentity(testCase)
for dScale = [1e-8, 0.3, -2.1]
    dVector = dScale*[1; -0.4; 0.2];
    [~, dJacobian] = RotationVectorToDCM(dVector);
    dExpected = integral(@(dTime) expm(dTime*skewSymm(dVector)), ...
        0,1,'ArrayValued',true,'AbsTol',1e-13);
    verifyEqual(testCase, dJacobian, dExpected, 'AbsTol', 2e-14);
end
end
