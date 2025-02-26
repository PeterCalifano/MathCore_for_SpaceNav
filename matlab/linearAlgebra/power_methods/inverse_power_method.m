function [v, lambda] = inverse_power_method(K, M, tol) %#codegen
%% PROTOTYPE
% [v, lambda] = inverse_power_method(K, M, tol)
% -------------------------------------------------------------------------
%% DESCRIPTION
% Inverse Power method to find the first eigenpair of the generalized
% eigenvalue problem (K-w^2*M)u = 0. If K is detected to be singular by
% check condition on its rank, then shifting is applied
% -------------------------------------------------------------------------
%% INPUT
% K: [NxN] Stiffness matrix 
% M: [NxN] Mass matrix
% tol: [scalar] defines the required accuracy of the solution
% Note: when the difference between lambda at k and k-1 becomes less than
%       tol, the iteration stops.
% -------------------------------------------------------------------------
%% OUTPUT
% v: [Nx1]: first eigenvector of the system
% lambda: [1] first eigenvalue of the system
% -------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano
% -------------------------------------------------------------------------
%% CHANGELOG
% v1: basic Inverse Power method with/without shifting coded and checked
%     against MATLAB function eig(), 25/08/2022
% -------------------------------------------------------------------------
%% Next upgrades

stop_condition = 1e3;
safe_var = 0;

% Dimension of the problem
[nrow, ncol] = size(K);

if rank(K) < min(nrow, ncol) % K is singular --> Shifting is applied
    sigma = 1000;

    K_tilde = K + sigma*M;

    % Generate a random vector as initial guess
    v0 = rand(nrow, 1);
    % Normalization to make the larget component of v_k equal to 1
    v0 = v0./max(v0);

    % First iteration
    v_hat0 = M*v0;
    lambda = 0;

    % Solve the linear system
    v = K_tilde\v_hat0;

    k = 1;
    disp(["Iteration number ", k]);
    err = 1;

    while err > tol && safe_var < stop_condition
        safe_var = safe_var + 1;

        k = k + 1;
        disp(["Iteration number ", k]); % for debug

        v_prev = v;
        lambda_prev = lambda;

        % First step
        v_prevhat = M*v_prev;

        % Second step: solve the linear system
        v = K\v_prevhat;

        lambda = (v'*K_tilde*v)/(v'*M*v); %(v_temp'*v_temp)/(v_temp'*v_prev);
        v = v./max(v);

        err = abs(lambda - lambda_prev);

    end

    lambda = lambda - sigma;

else

    % Generate a random vector as initial guess
    v0 = rand(nrow, 1);
    % Normalization to make the larget component of v_k equal to 1
    v0 = v0./max(v0);

    % First iteration
    v_hat0 = M*v0;
    lambda = 0;

    % Solve the linear system
    v = K\v_hat0;

    k = 1;
    disp(["Iteration number ", k]);
    err = 1;

    while err > tol
        safe_var = safe_var + 1;

        if safe_var == stop_condition
            warning('Iteration number limit reached. Solver has not reached convergence')
            break;
        end

        k = k + 1;
        disp(["Iteration number ", k]); % for debug

        v_prev = v;
        lambda_prev = lambda;

        % First step
        v_prevhat = M*v_prev;

        % Second step: solve the linear system
        v = K\v_prevhat;

        lambda = (v'*K*v)/(v'*M*v); %(v_temp'*v_temp)/(v_temp'*v_prev);
        v = v./max(v);

        err = abs(lambda - lambda_prev);

    end

    % Output of the while cycle: smallest eigenpair [v, lambda]

end

end