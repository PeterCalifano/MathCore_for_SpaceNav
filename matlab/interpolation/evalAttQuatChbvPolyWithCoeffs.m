function [dChbvInterpVector] = evalAttQuatChbvPolyWithCoeffs(ui32PolyDeg, ...
    ui32OutputSize, ...
    dEvalPoint, ...
    dChbvCoeffs, ...
    dSwitchIntervals, ...
    dDomainLB, ...
    dDomainUB, ...
    ui32PolyMaxDeg) %#codegen
arguments
    ui32PolyDeg         (1, 1) uint32   % {isscalar, isnumeric} % Commented for speed-up
    ui32OutputSize      (1, 1) uint32   % {isscalar, isnumeric}
    dEvalPoint          (1, 1) double   % {isscalar, isnumeric}
    dChbvCoeffs         (:, 1) double   % {isnumeric, ismatrix}
    dSwitchIntervals    (:, 2) double   % {isnumeric, ismatrix}
    dDomainLB           (1, 1) double   % {isscalar, isnumeric}
    dDomainUB           (1, 1) double   % {isscalar, isnumeric}
    ui32PolyMaxDeg      (1, 1) uint32   = ui32PolyDeg % {isscalar, isnumeric} 
end
%% PROTOTYPE
% [dChbvInterpVector] = evalAttQuatChbvPolyWithCoeffs(ui32PolyDeg, ...
%                                                       ui32OutputSize, ...
%                                                       dEvalPoint, ...
%                                                       dChbvCoeffs, ...
%                                                       dSwitchIntervals, ...
%                                                       dDomainLB, ...
%                                                       dDomainUB,
%                                                       ui32PolyMaxDeg)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg         (1, 1) uint32   % {isscalar, isnumeric} % Commented for speed-up
% ui32OutputSize      (1, 1) uint32   % {isscalar, isnumeric}
% dEvalPoint          (1, 1) double   % {isscalar, isnumeric}
% dChbvCoeffs         (:, 1) double   % {isnumeric, ismatrix}
% dSwitchIntervals    (:, 2) double   % {isnumeric, ismatrix}
% dDomainLB           (1, 1) double   % {isscalar, isnumeric}
% dDomainUB           (1, 1) double   % {isscalar, isnumeric}
% ui32PolyMaxDeg      (1, 1) uint32   = ui32PolyDeg % {isscalar, isnumeric}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-05-2024    Pietro Califano     First version, modified from general purpose utility. Validated.
% 18-07-2025    Pietro Califano     Fix basis and fitting problem errors
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% ui32PtrToLastCoeff = ui32PolyDeg * ui32OutputSize;
assert(length(dChbvCoeffs) == (ui32PolyMaxDeg+1)*ui32OutputSize, ...
    'Number of coefficients does not match output vector size.')

% Variables declaration
dChbvPolynomial     = zeros(ui32PolyMaxDeg + 1, 1);
dChbvInterpVector   = zeros(ui32OutputSize, 1);

% Compute scaled evaluation point
dScaledPoint = (2 * dEvalPoint - (dDomainLB + dDomainUB)) / (dDomainUB - dDomainLB);

% Get evaluated Chebyshev polynomials at scaled point
dChbvPolynomial(1:ui32PolyDeg+1) = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg);

% Compute interpolated output value by inner product with coefficients matrix
dChbvInterpVector(1:ui32OutputSize) = transpose( reshape(dChbvCoeffs,...
                                         ui32PolyDeg+1, ui32OutputSize) ) * dChbvPolynomial;

% DEVNOTE: not needed
% Switch sign of the interpolated value if required
% Check if within "switch intervals
% for idCheck = 1:size(dSwitchIntervals, 1) % TODO verify if this is allowed or iterable must be fixed
%     if dEvalPoint >= dSwitchIntervals(idCheck, 1) && dEvalPoint < dSwitchIntervals(idCheck, 2)
%         dChbvInterpVector(1:ui32OutputSize) = - dChbvInterpVector(1:ui32OutputSize);
%     end
% end

end
