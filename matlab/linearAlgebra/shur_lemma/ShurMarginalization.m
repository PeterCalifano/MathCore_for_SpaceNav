function [dMarginalCov] = ShurMarginalization(dCovMatrix, ...
                                              ui16FirstMargStateIdx) %#codegen
arguments
    dCovMatrix
    ui16FirstMargStateIdx
end
%% PROTOTYPE
% [dMarginalCov] = ShurMarginalization(dCovMatrix, ui16FirstMargStateIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computation of marginal covariance matrix from a Joint PDF covariance using Shur Complement. Correlations 
% are removed from the prior covariance of the marginalized states. The algorithm can work only on the
% bottom portion of the state vector.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dCovMatrix
% ui16FirstMargStateIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dMarginalCov
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 21-04-2024        Pietro Califano         First simple version coded.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
if coder.target('MATLAB') || coder.target('MEX')
    assert(size(dCovMatrix, 1) == size(dCovMatrix, 2));
    assert(ui16FirstMargStateIdx <= size(dCovMatrix, 1));
    assert(isscalar(ui16FirstMargStateIdx));
end

% Extract prior covariance of the states to marginalize, cross-correlations and remaining states covariance
LAMBDA22 = dCovMatrix(ui16FirstMargStateIdx:end, ui16FirstMargStateIdx:end);
LAMBDA12 = dCovMatrix(1:ui16FirstMargStateIdx-1, ui16FirstMargStateIdx:end);
LAMBDA11 = dCovMatrix(1:ui16FirstMargStateIdx-1, 1:ui16FirstMargStateIdx-1);

% Remove correlation terms from prior covariance to get marginalized covariance of the states
dMarginalCov = coder.nullcopy(zeros( size(dCovMatrix, 1) - ui16FirstMargStateIdx+1) );
dMarginalCov(:, :) = LAMBDA22 - transpose(LAMBDA12) * ( LAMBDA11\LAMBDA12 );


end
