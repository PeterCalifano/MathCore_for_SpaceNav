function [dRotatedVecArray] = Rot3dVecAboutDir(dRotDirArray, dVecToRotateArray, dRotAngleArray) %#codegen
arguments
    dRotDirArray            (3,:) double {ismatrix}
    dVecToRotateArray       (3,:) double {ismatrix}
    dRotAngleArray          (1,:) double {isvector}
end
%% PROTOTYPE
% [dRotatedVec] = Rot3dVecAboutDir(dRotDir, dVecToRotate, dRotAngle)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Rotation of a 3D vector dVecToRotate about a unit vector dRotDir by an angle dRotAngle, applying
% Rodriguez's formula. By construction, the formula works for any vector not parallel to the rotation.
% This is because it rotates the perpendicular component. NOTE: rotation is according to right-hand rule.
% REFERENCE:
% [1] J. Solà,  Quaternion kinematics for the error-state Kalman filter, Nov. 2017, page 19
% -------------------------------------------------------------------------
%% INPUT
% dRotDir:     [3x1]  Unit vector around which rotation takes place. The
%                       direction influences the sign of the rotation
%                       (i.e., the sign of this vector matters).
% dVec2Rotate: [3x1]  Vector to rotate about dRotDir
% dRotAngle:   [1]    [rad] Angle of which the component of dVec2Rotate
%                       parallel to dRotDir is rotated.
% -------------------------------------------------------------------------
%% OUTPUT
% dRotatedVec: [3x1]  Rotated vector
% -------------------------------------------------------------------------
%% CHANGELOG
% 17-10-2023    Pietro Califano    Retrieved from OM codes and re-worked.
% 16-03-2025    Pietro Califano    Upgrade to support vectorized operations
% -------------------------------------------------------------------------

%% Function code
ui32NumOfVectors = size(dVecToRotateArray, 2);

if coder.target("MATLAB") || coder.target("MEX")
    assert(ui32NumOfVectors == size(dRotDirArray, 2) || size(dRotDirArray, 2) == 1, ...
        'ERROR: invalid number of direction of rotation. Must be either 1 or equal to the number of vectors to rotate.')
    assert(ui32NumOfVectors == length(dRotAngleArray) || length(dRotAngleArray) == 1, ...
        'ERROR: invalid number of direction of rotation. Must be either 1 or equal to the number of vectors to rotate.')
end

% Handle special "conveniency" cases
if length(dRotAngleArray) == 1 % Explicit expansion of rotation angles
    dRotAngleArray = dRotAngleArray * ones(1, ui32NumOfVectors);
end

if size(dRotDirArray, 2) == 1 % Explicit expansion of rotation directions
    dRotDirArray = repmat(dRotDirArray, 1, ui32NumOfVectors);
end

% Perform rotations
dCosAngle = cos(dRotAngleArray); % Must be a row vector

% Rotate vector by keeping parallel component fixed and rotating perpendicular of dRotAngle
% (rotation is positive counterclockwise in the direction of dRotDir)
dRotatedVecArray = dVecToRotateArray .* dCosAngle + cross(dRotDirArray, dVecToRotateArray, 1) .* sin(dRotAngleArray) + ...
    dot(dRotDirArray, dVecToRotateArray, 1) .* (ones(1,length(dCosAngle)) - dCosAngle) .* dRotDirArray;

% Zero out numerical error
dRotatedVecArray(abs(dRotatedVecArray) < 2 * eps) = 0;
end
