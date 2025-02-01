function [i_dxKeplState_new] = keplerianProp(i_dxKeplState, i_dMu, i_dTstep, i_bConvert2Cart) %#codegen
%% PROTOTYPE
% [i_dxKeplState_new] = keplerianProp(i_dxKeplState, i_dMu, i_dTstep, i_bConvert2Cart)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% General purpose Keplerian propagator in Classical Keplerian elements. It
% propagates the input state vector forward of backward depending on
% i_dTstep input. The last element is interpreted as Eccentric or
% Hyperbolic anomaly depending on the eccentricity value. Check on SMA
% included, but negative sign expected if Hyperbolic.
% Implementation NOT robust in case of eccentricity near 1, only for
% Elliptical and Hyperbolic orbits (Circular and Parabolic NOT handled).
% REFERENCE:
% 1) Fundamentals of Astrodynamics and Applications - D. Vallado. Section
%    2.2, Kepler's Problem. Algorithms 3 and 4.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dxKeplState: [6, 1] Initial state vector in EA or HA
% i_dMu: [1] Gravitational parameter
% i_dTstep: [1] Time step of the propagation
% i_bConvert2Cart: [1] Flag for conversion to Cartesian
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% i_dxKeplState_new: [6, 1] Final state vector in EA or HA
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 19-10-2023    Pietro Califano     First prototype coded. VALIDATED.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Tunable tolerance and max iteration as input arguments
% 2) Conversion to Cartesian as optional execution
% -------------------------------------------------------------------------------------------------------------
%% Function code
EPS = 1e-6;
ECC_thr = 0.01; % Limit for check on eccentricity 

% Allocate output equal to input
i_dxKeplState_new = i_dxKeplState;

if not(exist('i_bConvert2Cart', 'var'))
    i_bConvert2Cart = false;
end


if abs(i_dTstep) > EPS
    % EXECUTE IF i_dTstep sufficiently large

    % Check state vector for execution path
    TOL = 1e-9; % [rad]
    ECC = i_dxKeplState(2);
    assert(ECC ~= 0, 'Circular orbit case not handled');
    MAX_ITER = 50;
    iterN = 0;

    if ECC > 1 + ECC_thr
        %% HYPERBOLIC ORBIT
        xHkep = i_dxKeplState;

        % Convert Hyperbolic anomaly to Mean anomaly
        % Mk = @(Hk, ECCk) ECCk *sinh(Hk) - Hk; Hk = xHkep(6).
        Mnext = ECC * sinh(xHkep(6)) - xHkep(6);

        % Update Mean anomaly for propagation
        Mnext = Mnext + sqrt(i_dMu/(-xHkep(1)^3)) * i_dTstep;

        ERR = 1; % Error for Newton-Raphson cycle

        % Initial guess generation
        if ECC < 1.6
            if (Mnext > -pi && Mnext < 0) || Mnext > pi
                Hk1 = Mnext - ECC;
            else
                Hk1 = Mnext + ECC;
            end
        else

            if ECC < 3.6 && abs(Mnext) > pi
                Hk1 = Mnext - sign(Mnext) * ECC;
            else
                Hk1 = Mnext/(ECC-1);
            end

        end

        % Newton-Raphson cycle
        while ERR > TOL && iterN <= MAX_ITER 

            % Execute iteration of Newton-Raphson to solve for H anomaly
            Hk2 = Hk1 +...
                ( Mnext - ECC*sinh(Hk1) + Hk1 )./(ECC * cosh(Hk1) - 1);

            % Evaluate error for convergence check
            ERR = abs(Hk2 - Hk1);
            Hk1 = Hk2;

            iterN = iterN + 1;

        end

        if iterN > MAX_ITER + 1
            warning(strcat("MAX_ITER = ", num2str(MAX_ITER) ," reached"));
        end

        % Assign output state
        i_dxKeplState_new(6) = Hk2;

    elseif ECC < 1 + ECC_thr && ECC > 0.0001
        %% ELLIPTICAL ORBIT
        xKep = i_dxKeplState;

        % Convert Eccentric anomaly to Mean anomaly
        % Mk = @(Ek, ECCk) Ek - ECC * sin(Ek) ; Ek = xKep(6).
        Mnext =  xKep(6) - ECC * sin(xKep(6));

        % Update Mean anomaly for propagation
        Mnext = Mnext + sqrt(i_dMu/(xKep(1)^3)) * i_dTstep;

        ERR = 1; % Error for Newton-Raphson cycle

        % Initial guess generation
        if (Mnext > -pi && Mnext < 0) || Mnext > pi
            Ek1 = Mnext - ECC;
        else
            Ek1 = Mnext + ECC;
        end

        % Newton-Raphson cycle
        while ERR > TOL && iterN <= MAX_ITER

            % Execute iteration of Newton-Raphson to solve for H anomaly
            Ek2 = Ek1 +...
                ( Mnext - Ek1 + ECC * sin(Ek1))./(1 - ECC*cos(Ek1));

            % Evaluate error for convergence check
            ERR = abs(Ek2 - Ek1);
            Ek1 = Ek2;

            iterN = iterN + 1;

        end

        if iterN >= MAX_ITER + 1
            warning(strcat("MAX_ITER = ", num2str(MAX_ITER) ," reached"));
        end

        % Assign output state
        i_dxKeplState_new(6) = Ek2;


    else
        % NEARLY PARABOLIC: NOT HANDLED
        assert(abs(ECC - 1) < ECC_thr, 'Eccentricity too close to 1')
    end

    if i_bConvert2Cart == true
        disp('Conversion NOT YET implemented')
    end

end
