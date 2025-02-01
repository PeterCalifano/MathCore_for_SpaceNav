function [tGrid, xState, info] = RKF45symCompatible(fDyn, x0, h0, t0, tf, order)
%% PROTOTYPE
% [tGrid, xState, info] = RK45(fDyn, x0, h0, t0, tf, order)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% 
% References:
% 1) Orbital Mechanics for Engineering Students Ed.4, Curtis
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% fDyn
% x0
% h0
% t0
% tf
% order
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% tGrid
% xState
% info
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-08-2023    Pietro Califano     RK4 integrator compatible with casadi
%                                   symbolic framework. Validated against
%                                   ode45 solution.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
% Future upgrade
%   1) Repackage function handle to take varargin?

% NvariableArgs = length(varargin);

% Check if symbolic objects are given as input
if not(isreal(x0))
    import casadi.* %#ok<SIMPT> 
    symbFlag = 1;
else
    symbFlag = 0;
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
C4 = [25.0/216.0, 0.0, 1408.8/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0];

% Coupling coefficients of RK5
C5 = [16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0];

% Difference of coupling coefficients for error control
% DeltaC = C5 - C4;
% DeltaC = - [-0.002777777777778; 0; 0.030253411306043; 0.029199893673578; -0.020000000000000; -0.036363636363636];

%% State vector info
Ndim = size(x0, 1);

%% Evaluate non-linear function stages
% Initialize arrays

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

odetimer = tic;
while tk < tf
    % BEGIN TIME STEP

    % Check Step size for final time
    if (tf - (tk + hk)) < 0
        % Adjusts last time step if it exceeds final time
        hk = tf - tk;
        notlastStep = false;
    end

    %% STAGES EVALUATION
    % Assign evaluation times of current timestep
    tkStages = tk + A.*hk;

    % Stage 1 evaluation
    fStages1 = fDyn(tkStages(1), xk);
    % Stage 2 evaluation
    fStages2 = fDyn(tkStages(2) , xk...
        + hk*B(2, 1)*fStages1);
    % Stage 3 evaluation
    fStages3 = fDyn(tkStages(3) , xk...
        + hk*B(3, 1)*fStages1...
        + hk*B(3, 2)*fStages2);
    % Stage 4 evaluation
    fStages4 = fDyn(tkStages(4) , xk...
        + hk*B(4, 1)*fStages1...
        + hk*B(4, 2)*fStages2...
        + hk*B(4, 3)*fStages3);
    % Stage 5 evaluation
    fStages5 = fDyn(tkStages(5) , xk...
        + hk*B(5, 1)*fStages1...
        + hk*B(5, 2)*fStages2...
        + hk*B(5, 3)*fStages3...
        + hk*B(5, 4)*fStages4);
    % Stage 6 evaluation
    fStages6 = fDyn(tkStages(6) , xk...
        + hk*B(6, 1)*fStages1...
        + hk*B(6, 2)*fStages2...
        + hk*B(6, 3)*fStages3...
        + hk*B(6, 4)*fStages4...
        + hk*B(6, 5)*fStages5);

    NfEval = NfEval + 6;


    if notlastStep
        % Pre-allocate storage for trajectory and time grid based on the
        % first hk that satisfies the tolerance (heuristic)
        tGridTemp = zeros(ceil((tf - t0)/h0), 1);
        
        if not(symbFlag)
            tGridTemp = zeros(ceil((tf-t0)/(h0)), 1);
            xStateTemp = zeros(Ndim, ceil((tf-t0)/(h0)) );
        end
    end


    % Update time
    tk = tk + hk;

    fStages = [fStages1, fStages2, fStages3, fStages4, fStages5, fStages6];

    % Evaluate next state at tk+1
    if order == 4
        xk = RK4eval(xk, fStages, hk, C4);
    elseif order == 5
        xk = RK5eval(xk, fStages, hk, C5);
    else
        error('Select either 4th or 5th order!')
    end

    % Update step counter
    stepNum = stepNum + 1;

    % Store time and trajectory
    tGridTemp(stepNum+1) = tk;

    if not(symbFlag)
        xStateTemp(:, stepNum+1) = xk;
    end
    % END TIME STEP
end

% Post-process output arrays
if stepNum > 1
    tGrid = [t0; tGridTemp(2:(stepNum+1), 1)];
else
    tGrid = [t0, tGridTemp(end)];
end

if symbFlag
    xState = xk;
else
    xState = [x0, xStateTemp(:, 2:(stepNum+1))];
end

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
%% RK4
    function xNext = RK4eval(xk, fStages, hk, C4)

        % Sum Stages
        fStagesSum = C4(1) * fStages(:, 1);

        for iC = 2:6
            fStagesSum = fStagesSum + C4(iC) * fStages(:, iC);
        end
        % Compute next value of xk
        xNext = xk + hk * fStagesSum;
    end
%% RK5
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