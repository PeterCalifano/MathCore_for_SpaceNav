function [o_dMarginalCov] = ShurMarginalization(i_dCovMatrix, i_ui16FirstMargStateIdx) %#codegen
%% PROTOTYPE
% [o_dMarginalCov] = ShurMarginalization(i_dCovMatrix, i_ui16FirstMargStateIdx) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computation of marginal covariance matrix from a Joint PDF covariance using Shur Complement. Correlations 
% are removed from the prior covariance of the marginalized states. The algorithm can work only on the
% bottom portion of the state vector.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dCovMatrix
% i_ui16FirstMargStateIdx
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dMarginalCov
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
assert(size(i_dCovMatrix, 1) == size(i_dCovMatrix, 2));
assert(i_ui16FirstMargStateIdx <= size(i_dCovMatrix, 1));
assert(isscalar(i_ui16FirstMargStateIdx));

% Extract prior covariance of the states to marginalize, cross-correlations and remaining states covariance
LAMBDA22 = i_dCovMatrix(i_ui16FirstMargStateIdx:end, i_ui16FirstMargStateIdx:end);
LAMBDA12 = i_dCovMatrix(1:i_ui16FirstMargStateIdx-1, i_ui16FirstMargStateIdx:end);
LAMBDA11 = i_dCovMatrix(1:i_ui16FirstMargStateIdx-1, 1:i_ui16FirstMargStateIdx-1);

% Remove correlation terms from prior covariance to get marginalized covariance of the states
o_dMarginalCov = coder.nullcopy(zeros( size(i_dCovMatrix, 1)-i_ui16FirstMargStateIdx+1) );
o_dMarginalCov(:, :) = LAMBDA22 - transpose(LAMBDA12) * ( LAMBDA11\LAMBDA12 );


end
