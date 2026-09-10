function dBasisMap = ComputeUnit3BasisTransport(dReference, dReferenceBasis, ...
    dDirection, dDirectionBasis) %#codegen
%% SIGNATURE
% dBasisMap = ComputeUnit3BasisTransport(dReference, dReferenceBasis, dDirection, dDirectionBasis)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Map tangent coordinates along the shortest sphere arc from a reference direction to a direction.
% Rotate the reference tangent vectors with the minimal rotation between the two directions, then
% express them in the destination basis. Both bases must be orthonormal and tangent to their inputs.
% Keeping the reference fixed lets callers change observation charts without path-dependent drift.
% Zero directions and antipodal directions are undefined and rejected by ComputeUnit3LocalError.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dReference         Reference direction; need not have unit length.
% dReferenceBasis    Two orthonormal tangent vectors at the reference direction.
% dDirection         Destination direction; need not have unit length.
% dDirectionBasis    Two orthonormal tangent vectors at the destination direction.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dBasisMap          Orthogonal 2-by-2 map from reference to destination tangent coordinates.
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 09-09-2026  Pietro Califano, Codex gpt-6    Add reference-based tangent-coordinate transport.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% ComputeUnit3LocalError, RotationVectorToDCM.
% -------------------------------------------------------------------------------------------------------------
arguments (Input)
    dReference      (3,1) double {mustBeFinite}
    dReferenceBasis (3,2) double {mustBeFinite}
    dDirection      (3,1) double {mustBeFinite}
    dDirectionBasis (3,2) double {mustBeFinite}
end
arguments (Output)
    dBasisMap       (2,2) double
end

% Use the sphere logarithm's small-angle handling and singularity checks.
[dLocalError, ~, dLogBasis] = ComputeUnit3LocalError(dReference, dDirection);
dReferenceUnit = dReference / norm(dReference);
dDirectionUnit = dDirection / norm(dDirection);

assert(norm(dReferenceBasis'*dReferenceBasis-eye(2),'fro') < 1e-10 && ...
    norm(dDirectionBasis'*dDirectionBasis-eye(2),'fro') < 1e-10 && ...
    norm(dReferenceBasis'*dReferenceUnit) < 1e-10 && ...
    norm(dDirectionBasis'*dDirectionUnit) < 1e-10, ...
    'ComputeUnit3BasisTransport:InvalidBasis', 'Bases must be orthonormal and tangent.');

% The logarithm is tangent to the reference; its normal cross product is the rotation vector.
dRotationVector = cross(dReferenceUnit, dLogBasis*dLocalError);
dRotation = RotationVectorToDCM(dRotationVector);
dBasisMap = dDirectionBasis' * dRotation * dReferenceBasis;

end
