function [dDCM, dLeftJacobian] = RotationVectorToDCM(dRotationVector, dSmallAngleThreshold) %#codegen
%% SIGNATURE
% [dDCM, dLeftJacobian] = RotationVectorToDCM(dRotationVector, dSmallAngleThreshold)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Evaluate R(phi) = Exp(skew(phi)) and optionally its SO(3) left Jacobian.
% The vector axis and norm define the rotation axis and angle. The differential
% satisfies R(phi+dphi) = Exp(skew(J_l(phi)*dphi))*R(phi) to first order.
% Rotation and Jacobian share the Rodrigues coefficients and their small-angle
% series. For passive rotations Exp(-skew(b)), evaluate this function at -b.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dRotationVector       Three-entry rotation vector [rad].
% dSmallAngleThreshold  Constant angle threshold for series evaluation [rad],
%                      default 1e-4. Existing one-output calls remain supported.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDCM           Rotation matrix associated with dRotationVector.
% dLeftJacobian  Additive-vector increment to left local rotation increment.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-06-2026    Pietro Califano     Add exact DCM construction from rotation vector.
% 09-09-2026  Pietro Califano, Codex gpt-6    Share coefficients with optional left Jacobian.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% skewSymm
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dRotationVector       (:,1) {mustBeNumeric}
    dSmallAngleThreshold  (1,1) {mustBeNumeric, mustBeNonnegative, coder.mustBeConst} = 1e-4
end
arguments (Output)
    dDCM          (3,3)
    dLeftJacobian (3,3)
end

%% Function code

assert(numel(dRotationVector) == 3, ...
    'RotationVectorToDCM:InvalidInput', ...
    'Rotation vector must contain exactly 3 elements.');

% Compiler constant
bComputeJacobian = coder.const(nargout > 1);

% Compute rotation angle and vector norm
dRotationVector = dRotationVector(:);
dRotationAngle = norm(dRotationVector);

% The limiting rotation and differential are both identity.
dLeftJacobian = eye(3);
if dRotationAngle <= eps
    dDCM = eye(3);
    return
end

% Compute rotation angle and vector skew matrix
dRotationAngleSq = dRotationAngle * dRotationAngle;
dSkewRotationVector = skewSymm(dRotationVector);
dSkewSquared = dSkewRotationVector * dSkewRotationVector;

if dRotationAngle <= eps + coder.const(dSmallAngleThreshold)
    % Series avoid subtraction of nearly equal trigonometric terms near zero.
    dSinCoefficient = 1.0 - dRotationAngleSq/6.0 + dRotationAngleSq^2/120.0;
    dCosCoefficient = 0.5 - dRotationAngleSq/24.0 + dRotationAngleSq^2/720.0;
    
    if bComputeJacobian
        dJacCoefficient = 1/6 - dRotationAngleSq/120 + dRotationAngleSq^2/5040;
    end
else
    % Full Rodrigues formula and left Jacobian
    dSinCoefficient = sin(dRotationAngle)/dRotationAngle;
    dCosCoefficient = (1.0-cos(dRotationAngle))/dRotationAngleSq;

    if bComputeJacobian
        dJacCoefficient = (1.0-dSinCoefficient)/dRotationAngleSq;
    end
end

dDCM = eye(3) + dSinCoefficient*dSkewRotationVector + dCosCoefficient*dSkewSquared;

if bComputeJacobian
    dLeftJacobian = eye(3) + dCosCoefficient*dSkewRotationVector + dJacCoefficient*dSkewSquared;
end
end
