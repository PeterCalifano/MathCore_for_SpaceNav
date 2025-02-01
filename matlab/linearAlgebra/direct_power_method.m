function [v, lambda] = direct_power_method(A, tol) %#codegen
%% PROTOTYPE
% [v, lambda] = direct_power_method(A, tol)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Direct Power method to find the dominant eigenpair of matrix A
% -------------------------------------------------------------------------
%% INPUT
% A: [NxN] matrix describing the eigenvalue problem
% tol: [scalar] defines the required accuracy of the solution
% Note: when the difference between lambda at k and k-1 becomes less than
%       tol, the process stops
% -------------------------------------------------------------------------
%% OUTPUT
% u: [Nx1]: dominant eigenvector of the system
% lambda: [1] dominant eigenvalue of the system
% -------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano
% -------------------------------------------------------------------------
%% CHANGELOG
% V1: method implementation as standalone function, 05/12/2021
% -------------------------------------------------------------------------
%% Next upgrades
%


%% Function code
stop_condition = 1e8;
safe_var = 0;

% Dimension of the problem
[nrow, ~] = size(A);

% Generate a random vector as initial guess
v0 = rand(nrow, 1);
v0 = v0./norm(v0);

% First iteration
v = A*v0;
lambda = 0;

k = 1;
err = 1;

while err > tol && safe_var < stop_condition
    safe_var = safe_var + 1; 

    disp(["Iteration number ", k]); % for debug

    v_prev = v;
    lambda_prev = lambda;

    v_temp = A*v; 
    lambda = (v_temp'*v_temp)/(v_temp'*v_prev);
    v = v_temp./norm(v_temp);

    err = abs(lambda - lambda_prev);
    k = k + 1;
end


end