function [o_dRotatedVec] = Rot3dVecAboutDir(i_dRotDir, i_dVecToRotate, i_dRotAngle) %#codegen
%% PROTOTYPE
% [o_dRotatedVec] = Rot3dVecAboutDir(i_dRotDir, i_dVecToRotate, i_dRotAngle)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Rotation of a 3D vector i_dVecToRotate about a unit vector i_dRotDir 
% by an angle i_dRotAngle, applying Rodriguez's formula.
% NOTE: rotation is according to right-hand rule.
% REFERENCE:
% [1] J. Solà,  Quaternion kinematics for the error-state Kalman filterp, Nov. 2017, page 19
% -------------------------------------------------------------------------
%% INPUT
% i_dRotDir:     [3x1]  Unit vector around which rotation takes place. The 
%                       direction influences the sign of the rotation 
%                       (i.e., the sign of this vector matters).
% i_dVec2Rotate: [3x1]  Vector to rotate about i_dRotDir 
% i_dRotAngle:   [1]    [rad] Angle of which the component of i_dVec2Rotate
%                       parallel to i_dRotDir is rotated.
% -------------------------------------------------------------------------
%% OUTPUT
% o_dRotatedVec: [3x1]  Rotated vector
% -------------------------------------------------------------------------
%% CHANGELOG
% 17-10-2023    Pietro Califano    Retrieved from OM codes and re-worked.
% -------------------------------------------------------------------------


%% Function code
cosAngle = cos(i_dRotAngle);
o_dRotatedVec = i_dVecToRotate * cosAngle + cross(i_dRotDir, i_dVecToRotate) * sin(i_dRotAngle) + ...
            i_dRotDir * dot(i_dRotDir, i_dVecToRotate) * (1-cosAngle);

end
