function o_dLnew = cholRank1Update(i_dL, i_dU, i_dBeta, i_dAlpha) %#codegen
arguments
i_dL     (:, :) double {ismatrix, isnumeric}
i_dU     (:, :) double {ismatrix, isnumeric}
i_dBeta  (1, 1) double {isnumeric, isscalar} = 1.0
i_dAlpha (1, 1) double {isnumeric, isscalar} = 1.0
end
%% PROTOTYPE
% o_dLnew = cholRank1Update(i_dL, i_dU, i_dBeta, i_dAlpha)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing the Cholesky Rank 1 Update on the columns of matrix
% i_dU using "Triangular Rank-One update". Given a matrix Sigma with
% Cholesky (lower) factor L, a vector v and two real scalar alpha and beta,
% the algorithm computes the updated Cholesky factor Lnew such that:
% Lnew = chol(SigmaNew), where SigmaNew = alpha*Sigma + beta * v*v';
% REFERENCE:
% 1) A More Efficient Rank-one Covariance Matrix Update for Evolution 
% Strategies, Krause, Igel, 2015
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dL     (:, :)               Original Cholesky factor L (can be either lower or upper, automatically handled)
% i_dU     (:, :)               Update Matrix with update vectors v as columns
% i_dBeta  (1, 1) double = 1.0  Scaling coefficient of Update matrix
% i_dAlpha (1, 1) double = 1.0  Scaling coefficient of cholesky factor
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dLnew
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 17-09-2023    Pietro Califano     First prototype coded. Debug needed:
%                                   slightly output with respect to chol().
%                                   Error is in the order of 10^-2 on the
%                                   diagonal whereas other elements are
%                                   below machine precision in double.
% 04-02-2024    Pietro Califano     Code validated against MATLAB cholUpdate
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Check if input Cholesky factor is lower triangular, else transpose
UPPER = logical(istriu(i_dL));

if UPPER
    dL = transpose(i_dL);
else
    dL = i_dL;
end

% Check if scaling is present
if nargin < 4
    i_dAlpha = 1.0;
    if nargin < 3
        i_dBeta = 1.0;
    end
end

% Determine number of columns of U
NrowColU = uint16(size(i_dU));
NrowU = NrowColU(1);
NcolU = NrowColU(2);

N = uint16(size(dL, 1));
assert(NrowU == N, 'Dimension of rows of L and U does not match!')

% Allocate output matrix
dLnew = zeros(size(dL));

% Input: Cholesky factor i_dL of size N
alphaSqrt = sqrt(i_dAlpha);

for idCol = 1:NcolU

    % Triangular Rank1 Update
    v = i_dU(:, idCol);
    w = v; % Initialize rank-1 vector
    b = 1.0; % Initialize b coefficient

    for idj = 1:N
        % Compute entry of new Cholesky factor
        dLnew(idj, idj) = sqrt( i_dAlpha * dL(idj, idj)^2 ...
            + (i_dBeta/b) * w(idj)^2 );
        
        gamma = i_dAlpha * dL(idj, idj)^2 * b + i_dBeta * w(idj)^2;

        for idk = idj+1:N     
            % Update kth value in "update" vector
            w(idk) = w(idk) - (w(idj)/dL(idj, idj)) * alphaSqrt * dL(idk, idj);
            % Compute entry of new Cholesky factor
            dLnew(idk, idj) = alphaSqrt * (dLnew(idj, idj)/dL(idj, idj)) * dL(idk, idj) +...
                (dLnew(idj, idj) * i_dBeta * w(idj)) * w(idk)/gamma;         
        end

        % Update auxiliary variable b
        b = b + i_dBeta * w(idj)^2/(i_dAlpha * dL(idj, idj)^2);

    end
end

% If input was Upper, provide output as Upper Cholesky factor
if UPPER
    o_dLnew = transpose(dLnew);
else
    o_dLnew = dLnew;
end
