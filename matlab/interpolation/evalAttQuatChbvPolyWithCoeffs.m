function [dChbvInterpVector] = evalAttQuatChbvPolyWithCoeffs(ui32PolyDeg, ui32OutputSize, ...
    dEvalPoint, dChbvCoeffs, dswitchIntervals, dDomainLB, dDomainUB) %#codegen
arguments
    ui32PolyDeg         (1, 1) uint32   {isscalar, isnumeric}
    ui32OutputSize      (1, 1) uint32   {isscalar, isnumeric}
    dEvalPoint          (1, 1) double   {isscalar, isnumeric}
    dChbvCoeffs         (:, 1) double
    dswitchIntervals    (:, 2) double 
    dDomainLB           (1, 1) double   {isscalar, isnumeric}
    dDomainUB           (1, 1) double   {isscalar, isnumeric}
end
% TODO (PC) UPDATE IMPLEMENTATION
%% PROTOTYPE
% [o_dChbvInterpVector] = evalChbvPolyWithCoeffs(i_ui8PolyDeg, i_ui8OutputSize, ...
    % i_dEvalPoint, i_dChbvCoeffs, i_dswitchIntervals, i_dDomainLB, i_dDomainUB)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg
% i_ui8OutputSize
% i_dEvalPoint
% i_dChbvCoeffs
% i_bIsSignSwitched
% i_dDomainLB
% i_dDomainUB
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-05-2024        Pietro Califano         First version, modified from general purpose utility. Validated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

assert(length(dChbvCoeffs) == ui32PolyDeg*ui32OutputSize, ...
    'Number of coefficients does not match output vector size.')

% Variables declaration
dChbvPolynomial = coder.nullcopy(zeros(ui32PolyDeg+1, 1));
dChbvInterpVector = coder.nullcopy(zeros(ui32OutputSize, 1));

% Compute scaled evaluation point
dScaledPoint = (2 * dEvalPoint - (dDomainLB+dDomainUB)) / (dDomainUB-dDomainLB);
% Get evaluated Chebyshev polynomials at scaled point
dChbvPolynomial(:) = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint);

% Compute interpolated output value by inner product with coefficients matrix
dChbvInterpVector(:) = transpose( reshape(dChbvCoeffs,...
    ui32PolyDeg, ui32OutputSize) ) * dChbvPolynomial(2:end);

% Switch sign of the interpolated value if required
% Check if within "switch intervals
for idCheck = 1:size(dswitchIntervals, 1)
    if dEvalPoint >= dswitchIntervals(idCheck, 1) && dEvalPoint < dswitchIntervals(idCheck, 2)
        dChbvInterpVector(:) = - dChbvInterpVector(:);
    end
end

end
