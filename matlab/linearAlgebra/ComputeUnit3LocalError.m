function [dLocalVec2, dOtherJacobian, dTangentBasis] = ...
    ComputeUnit3LocalError(dUnit3Base, dUnit3Other) %#codegen
%% SIGNATURE
% [dLocalVec2, dOtherJacobian, dTangentBasis] = ComputeUnit3LocalError(dUnit3Base, dUnit3Other)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Express the sphere logarithm in a deterministic orthonormal basis at the base direction.
% Differentiate with respect to the unnormalized other vector, holding the base and basis fixed.
% Zero vectors and directions within the antipodal guard have no supported chart and are rejected.
% Implemented following GTSAM Unit3::localCoordinates implementation.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dUnit3Base       Nonzero base vector; its magnitude is ignored.
% dUnit3Other      Nonzero vector to map into the base tangent plane.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dLocalVec2       Two angular coordinates [rad].
% dOtherJacobian   Derivative with respect to dUnit3Other, including normalization.
% dTangentBasis    Two orthonormal columns perpendicular to the normalized base.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Add analytic derivatives and reject undefined charts.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% Build3dOrthonormalBasis.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dUnit3Base  (3,1) double {mustBeFinite}
    dUnit3Other (3,1) double {mustBeFinite}
end
arguments (Output)
    dLocalVec2     (2,1) double
    dOtherJacobian (2,3) double
    dTangentBasis  (3,2) double
end

% Normalization removes the radial degree of freedom before taking the logarithm.
dBaseNorm = norm(dUnit3Base);
dOtherNorm = norm(dUnit3Other);

assert(dBaseNorm > 0 && dOtherNorm > 0, 'ComputeUnit3LocalError:ZeroVector', ...
    'Sphere coordinates require two nonzero vectors.');

dBase = dUnit3Base/dBaseNorm;
dOther = dUnit3Other/dOtherNorm;

[dAxis1, dAxis2] = Build3dOrthonormalBasis(dBase);
dTangentBasis = [dAxis1, dAxis2];
dCosAngle = min(1.0, max(-1.0, dot(dBase,dOther)));

assert(dCosAngle > -1.0 + 1e-12, 'ComputeUnit3LocalError:Antipodal', ...
    'The sphere logarithm is undefined near the antipodal direction.');

% Use the first order approximation series near zero angle to avoid cancellation in the scale derivative.
dCosGap = 1.0 - dCosAngle;

if dCosGap < 1e-4
    dScale = 1.0 + dCosGap*(1/3 + dCosGap*(2/15 + dCosGap*2/35));
    dScaleDerivative = -1/3 - dCosGap*(4/15 + dCosGap*6/35);
else
    dSinAngle = norm(cross(dBase, dOther));
    dScale = atan2(dSinAngle, dCosAngle) / dSinAngle;
    dScaleDerivative = (dCosAngle*dScale-1.0) / (dSinAngle*dSinAngle);
end

dTangentVector = dTangentBasis' * (dOther-dCosAngle*dBase);
dLocalVec2 = dScale * dTangentVector;
dOtherJacobian = zeros(2,3);

if coder.const(nargout > 1)
    dNormalizeJacobian = ( eye(3) - dOther*dOther' ) / dOtherNorm;
    dOtherJacobian = ( dScale * dTangentBasis' + ...
        dScaleDerivative * dTangentVector * dBase' ) * dNormalizeJacobian;
end
end
