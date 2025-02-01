function [dChbvCoeffs, dScaledInterpDomain, strfitStats] = fitChbvPolynomials(ui32PolyDeg, ...
    dInterpDomain, ...
    dDataMatrix, ...
    dDomainLB, ...
    dDomainUB, ...
    bENABLE_FIT_CHECK) %#codegen
arguments
    ui32PolyDeg         (1, 1) uint32
    dInterpDomain       (:, 1) double
    dDataMatrix         (:, :) double
    dDomainLB           (1, 1) double 
    dDomainUB           (1, 1) double  
    bENABLE_FIT_CHECK   (1, 1) logical = true
end
%% PROTOTYPE
% [dChbvCoeffs, dScaledInterpDomain, strfitStats] = fitChbvPolynomials(ui32PolyDeg, ...
%     dInterpDomain, ...
%     dDataMatrix, ...
%     dDomainLB, ...
%     dDomainUB, ...
%     bENABLE_FIT_CHECK) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function for fitting interpolation coefficients of Chebyshev Polynomial up to the specified degree. The
% interpolant maps the input 1D domain to a N-dimensional domain as specified by the 1st dimension of the
% data matrix. However, note that each entry of the jth sample is interpolated by a different polynomial.
% The output is a 1D vector containing the coefficients for Chebyshev Polynomials from the 1st to PolyDeg-th
% degree. The scaling to [-1,1] domain is automatically handled.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg
% i_dInterpDomain
% i_dDataMatrix
% i_dDomainLB
% i_dDomainUB
% i_bENABLE_AUTO_CHECK
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dChbvCoeffs
% o_dScaledInterpDomain
% o_strfitStats
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-04-2024        Pietro Califano         First version. Validated.
% 08-05-2024        Pietro Califano         Updated with error checks.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

if nargin < 6
    bENABLE_FIT_CHECK = true;
end

% Check input dimensions
% i_dDataMatrix: [L, N] where N is the number of points, L is the output vector size
assert(size(dDataMatrix, 2) == length(dInterpDomain));
assert(ui32PolyDeg > 2);


assert(length(dInterpDomain) >= ui32PolyDeg +1);

% Get size of the output vector
ui8OutputSize = size(dDataMatrix, 1);

% Allocate output matrix
dChbvCoeffs = zeros(ui8OutputSize*(ui32PolyDeg), 1);


if nargin < 3
    dDomainUB = max(dInterpDomain, [], 'all');
    dDomainLB = min(dInterpDomain, [], 'all');
end

% Compute scaled domain
dScaledInterpDomain = (2.*dInterpDomain - (dDomainUB+dDomainLB))./(dDomainUB-dDomainLB);

% Compute regressors matrix on scaled domain
dRegrMatrix = zeros(ui32PolyDeg, size(dDataMatrix, 2));

for idN = 1:size(dDataMatrix, 2)

    % Evaluate Chebyshev polynomial at scaled point
    tmpChbvPoly = EvalRecursiveChbv(ui32PolyDeg, dScaledInterpDomain(idN));
    dRegrMatrix(:, idN) = tmpChbvPoly(2:end);
    
end

% Compute fit coefficients matrix (transposed)
% Xmat = Cmat * Phi: [LxN] = [LxM]*[MxN] where M: poly degree, N: number of samples, L: output vector size
% ith Chebyshev polynomial i=1,...N, along each column
% jth element of ith sample has coefficients along each row of Cmat
dChbvCoeffs_matrixT = dRegrMatrix' \ dDataMatrix'; % Solve the transposed problem
% Flatten matrix to 1D vector
dChbvCoeffs(1:end) = dChbvCoeffs_matrixT(:);

if bENABLE_FIT_CHECK == true && not(isempty('evalChbvPolyWithCoeffs.m'))
        [strfitStats] = checkFitChbvPoly(ui32PolyDeg, dInterpDomain, dChbvCoeffs, ...
            dDataMatrix, dDomainLB, dDomainUB, false);
else
    strfitStats = struct();
end

end
