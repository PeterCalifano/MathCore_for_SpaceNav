function [dChbvInterpVector] = evalChbvPolyWithCoeffs(ui32PolyDeg, ...
                                                      ui32OutputSize, ...
                                                      dEvalPoint, ...
                                                      dChbvCoeffs, ...
                                                      dDomainLB, ...
                                                      dDomainUB, ...
                                                      ui32PtrToLastCoeff, ...
                                                      ui32PolyMaxDeg) %#codegen
arguments
    ui32PolyDeg         (1,1) uint32    % {isscalar, isnumeric} % Commented for speed-up
    ui32OutputSize      (1,1) uint32    % {isscalar, isnumeric}
    dEvalPoint          (1,1) double    % {isscalar, isnumeric}
    dChbvCoeffs         (:,1) double    % {isnumeric, ismatrix}
    dDomainLB           (1,1) double    % {isscalar, isnumeric}
    dDomainUB           (1,1) double    % {isscalar, isnumeric}
    ui32PtrToLastCoeff  (1,1) uint32    % {isscalar, isnumeric} = length(dChbvCoeffs)
    ui32PolyMaxDeg      (1,1) uint32    % {isscalar, isnumeric} = ui32PolyDeg
end
%% PROTOTYPE
% TODO: update doc
% [dChbvInterpVector] = evalChbvPolyWithCoeffs(ui8PolyDeg, ...
%                                              ui8OutputSize, ...
%                                              dEvalPoint, ...
%                                              dChbvCoeffs, ...
%                                              dDomainLB, ...
%                                              dDomainUB) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% ui32PolyDeg         (1,1) uint32    % {isscalar, isnumeric} % Commented for speed-
% ui32OutputSize      (1,1) uint32    % {isscalar, isnumeric}
% dEvalPoint          (1,1) double    % {isscalar, isnumeric}
% dChbvCoeffs         (:,1) double    % {isnumeric, ismatrix}
% dDomainLB           (1,1) double    % {isscalar, isnumeric}
% dDomainUB           (1,1) double    % {isscalar, isnumeric}
% ui32PtrToLastCoeff  (1,1) uint32    % {isscalar, isnumeric} = length(dChbvCoeffs)
% ui32PolyMaxDeg      (1,1) uint32    % {isscalar, isnumeric} = ui32PolyDeg
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-04-2024    Pietro Califano     First version, verified in unit test.
% 01-02-2025    Pietro Califano     Upgrade of functions for codegen with static-sized arrays
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Check length of used coefficients 
assert(ui32PtrToLastCoeff == ui32PolyDeg * ui32OutputSize, ...
    'Number of coefficients does not match the expected size!')

% Variables definition (static size)
dChbvPolynomial     = coder.nullcopy(zeros(ui32PolyMaxDeg + 1, 1));
dChbvInterpVector   = coder.nullcopy(zeros(ui32OutputSize, 1));

% Compute scaled evaluation point 
dScaledPoint = coder.nullcopy(0.0);
dScaledPoint(:) = (2 * dEvalPoint - (dDomainLB + dDomainUB)) / (dDomainUB - dDomainLB); % scalar

% Get evaluated Chebyshev polynomials at scaled point
dChbvPolynomial(1:ui32PolyDeg+1) = EvalRecursiveChbv(ui32PolyDeg, dScaledPoint, ui32PolyMaxDeg); % TODO 

% Compute interpolated output value by inner product with coefficients matrix
dChbvInterpVector(1:ui32OutputSize) = transpose( reshape(dChbvCoeffs(1:ui32PtrToLastCoeff), ui32PolyDeg, ui32OutputSize) ) * ...
                                                    dChbvPolynomial(2:ui32PolyDeg + 1);


end
