function [tGrid, xState, info] = PropagateRKF45v2(fDyn, x0, h0, t0, tf, tol) %#codegen
%% PROTOTYPE
% [tGrid, xState, info] = RKF45(fDyn, x0, h0, t0, tf, tol)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Integrator function implementing Runge-Kutta-Fehlberg embedded formulas 
% of order 4-5 with adaptive law to determine the step size, by comparing
% the solutions of the two orders. The function is able to provide the
% entire trajectory as computed by the solver, but an initial memory
% overhead is necessary for pre-allocation. Interpolation not supported.
% References:
% 1) Orbital Mechanics for Engineering Students Ed.4, Curtis
% 2)
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% fDyn, x0, h0, t0, tf, tol
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tGrid, xState, info
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 30-07-2023    Pietro Califano     Coefficients hardcoded for RK45
% 31-07-2023    Pietro Califano     Coding of integrator completed.
%                                   Validation and interface for generic
%                                   RHS still not completed.
% 06-08-2023    Pietro Califano     Interface for input functions created.
%                                   First release version.
% 08-08-2023    Pietro Califano     Upgrade of RK5 eval with more clever
%                                   initialization of variables.
% 04-09-2023    Pietro Califano     Modified to make operational flow more
%                                   similar to reference implementation.
%                                   Accuracy is comparable to ode45, but
%                                   speed is much lower.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
% Future upgrade
%   1) Repackage function handle to take varargin?

% NvariableArgs = length(varargin);

if exist('tol', 'var')
    tol = 1e-8;
end


%% Function code
% Coefficients
% Weights of the solution at k time
A = [0.0; 1.0/4.0; 3.0/8.0; 12.0/13.0; 1.0; 1.0/2.0]; 

% Weights of the non linear function stages
B = [0.0,            0.0,           0.0,            0.0,            0.0;
     1.0/4.0,        0.0,           0.0,            0.0,            0.0;
     3.0/32.0,       9.0/32.0,      0.0,            0.0,            0.0;
     1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0,  0.0,            0.0; 
     439.0/216.0,   -8.0,           3680.0/513.0,   -845.0/4104.0,  0.0;
     -8.0/27.0       2.0,           -3544.0/2565.0  1859.0/4104.0,  -11.0/40.0];

% Coupling coefficients of RK4
% C4 = [25.0/216.0, 0.0, 1408.8/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0];

% Coupling coefficients of RK5
C5 = [16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0];

% Difference of coupling coefficients for error control
% DeltaC = C5 - C4;
DeltaC = - [-0.002777777777778; 0; 0.030253411306043; 0.029199893673578; -0.020000000000000; -0.036363636363636];

%% State vector info
Ndim = size(x0, 1);

%% Evaluate non-linear function stages
% Initialize arrays

fStages = zeros(Ndim, 6);

stepNum = 0;
NfEval = 0;

if not(exist('params', 'var'))
    params = []; %#ok<NASGU>
end

% Initialize integration
xk = x0;
tk = t0;
hk = h0;
notlastStep = true;
deltaStep = 1;

odetimer = tic;
while tk < tf
    % BEGIN TIME STEP

    % Check Step size for final time
    if (tf - (tk + deltaStep*hk)) < 0
        % Adjusts last time step if it exceeds final time
        hk = tf - tk;
        deltaStep = 1;
        notlastStep = false;
    end

    %% STAGES EVALUATION
    % Assign evaluation times of current timestep
    tkStages = tk + A.*deltaStep*hk;

    % Stage 1 evaluation
    fStages(:, 1) = fDyn(tkStages(1), xk);
    % Stage 2 evaluation
    fStages(:, 2) = fDyn(tkStages(2) , xk...
        + deltaStep*hk*B(2, 1)*fStages(:, 1));
    % Stage 3 evaluation
    fStages(:, 3) = fDyn(tkStages(3) , xk...
        + deltaStep*hk*B(3, 1)*fStages(:, 1)...
        + deltaStep*hk*B(3, 2)*fStages(:, 2));
    % Stage 4 evaluation
    fStages(:, 4) = fDyn(tkStages(4) , xk...
        + deltaStep*hk*B(4, 1)*fStages(:, 1)...
        + deltaStep*hk*B(4, 2)*fStages(:, 2)...
        + deltaStep*hk*B(4, 3)*fStages(:, 3));
    % Stage 5 evaluation
    fStages(:, 5) = fDyn(tkStages(5) , xk...
        + deltaStep*hk*B(5, 1)*fStages(:, 1)...
        + deltaStep*hk*B(5, 2)*fStages(:, 2)...
        + deltaStep*hk*B(5, 3)*fStages(:, 3)...
        + deltaStep*hk*B(5, 4)*fStages(:, 4));
    % Stage 6 evaluation
    fStages(:, 6) = fDyn(tkStages(6) , xk...
        + deltaStep*hk*B(6, 1)*fStages(:, 1)...
        + deltaStep*hk*B(6, 2)*fStages(:, 2)...
        + deltaStep*hk*B(6, 3)*fStages(:, 3)...
        + deltaStep*hk*B(6, 4)*fStages(:, 4)...
        + deltaStep*hk*B(6, 5)*fStages(:, 5));

    NfEval = NfEval + 6;

    % Evaluate solution with current step (delta * hk)
    xkErrEval = RK5eval(xk, fStages, deltaStep * hk, C5);

    %% ERROR EVALUATION
    % Compute error
    fERR = zeros(Ndim, 1);
    for iE = 1:6
        fERR = fERR + DeltaC(iE) * fStages(:, iE);
    end
    % Compute state error
    errVec = (deltaStep * hk) * fERR;
    % Compute maximum error with current (deltaStep * hk)
    maxErr = max(abs(errVec));
    % Compute the maximum allowed (relative)
    xMax = max(abs(xkErrEval)); % FAULT IS HERE: this must change to compute whether the relative tolerance is satisfied
    maxErr_allowed = tol * max(xMax, 1.0);

    %% ADAPTIVE STEP EVALUATION
    % From reference (2)

    if notlastStep
        % Skip error control if last step
        %         hNext = 0.9 * hk * (tol/maxErr)^(1/5);

        % Update timestep size
        %         hk = hNext;


        if maxErr > maxErr_allowed
            % NOTE: If error exceed relative tol:
            % Repeat error evaluation cycle and update deltaStep

            % Set repeat flag to update deltaStep
            repeatFlag = true;

            % Determine minimum step for satefy
            hmin = 16*eps(tk);
            deltaStep = (maxErr_allowed/(maxErr + eps)) ^ (1.0/5.0);


            if hk < hmin
                warning('Minimum step size reached. Tolerance can not be satisfied: minimum enforced.')
                hk = hmin;
                repeatFlag = false;
            end

        else 

            % Error tolerance satisfied, go ahead
            repeatFlag = false;
            if (tk - t0) < eps
                % Pre-allocate storage for trajectory and time grid based on the
                % first hk that satisfies the tolerance (heuristic)
                ALLOC_SIZE = ceil((tf-tk)/hk);

                tGridTemp = zeros(ALLOC_SIZE, 1);
                xStateTemp = zeros(Ndim, ALLOC_SIZE);
            end

        end

        % Modify step size for next iteration
        hk = min(deltaStep * hk, 4*hk);
        deltaStep = 1; % Reset Delta Step

    end


    %% MOVE TO NEXT STEP
    if repeatFlag == false

        % Update solution time
        tk = tk + hk;
        xk = xkErrEval;

        % Update step counter
        stepNum = stepNum + 1;

        % Store time and trajectory
        tGridTemp(stepNum+1) = tk;
        xStateTemp(:, stepNum+1) = xk;


    end

    % END TIME STEP
end

% Post-process output arrays 
tGrid = [t0; tGridTemp(2:(stepNum+1), 1)];
xState = [x0, xStateTemp(:, 2:(stepNum+1))];

elapsedTime = toc(odetimer);

% Store info if required
if nargout > 2
    info.stepNum = stepNum;
    info.NfEval = NfEval;
    info.minStep = min(abs(diff(tGrid)));
    info.maxStep = max(abs(diff(tGrid)));
    info.avgStep = mean(abs(diff(tGrid)));
    info.runtime = elapsedTime;
    info.medStep = median(abs(diff(tGrid)));
end

%% LOCAL FUNCTIONS

%     function xNext = RK4eval(xk, fStages, hk, C4)
%         xNext = xk + hk * C4 .* fStages;
%     end
%
    function xNext = RK5eval(xk, fStages, hk, C5)

        % Sum Stages
        fStagesSum = C5(1) * fStages(:, 1);

        for iC = 2:6
            fStagesSum = fStagesSum + C5(iC) * fStages(:, iC);
        end
        % Compute next value of xk
        xNext = xk + hk * fStagesSum;
    end

end



