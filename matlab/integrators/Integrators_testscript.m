close all
clear
clc

%% Integrators test script
% BASED ON EXERCISE 3 OF ASSIGNMENT 1

% Equation: xdot = -(t^2 - 1)*x/(t^2 + 1)
x0 = 1; % at time t0 = 0
x_analytical = @(t) exp(2*atan(t) - t);

% Time span bounds
t0 = 0;
tf = 2;

% RHS of the system
f = @(x, t) -(t^2 - 1)*x./(t^2 + 1);
% Step size
h = [0.5, 0.2, 0.05, 0.01];

theta = 0.4;
evalcounter_RK2 = nan(1, length(h));
evalcounter_RK4 = nan(1, length(h));

for idh = 4

    % Select h
    tspan = t0:h(idh):tf;

    tic
    [x_RK4, evalcounter_RK4(idh)] = RKN_integrator(f, tspan, x0, 4);
    toc

    tic
    [~, x_BI2_04, evalcounter_BI2_04] = BIn_theta_integrator(f, tspan, x0, 2, h(idh), theta);
    toc

    tic
    [~, x_AB3, evalcounter_AB3] = Multistep_integrator(f, tspan, x0, 'ab', 3);
    toc

    tic
    [~, x_AM3, evalcounter_AM3] = Multistep_integrator(f, tspan, x0, 'am', 3);
    toc

    tic
    [~, x_ABM3, evalcounter_ABM3] = Multistep_integrator(f, tspan, x0, 'abm', 3);
    toc

    tic
    [~, x_BDF3, evalcounter_BDF3] = Multistep_integrator(f, tspan, x0, 'bdf', 3);
    toc
end


% PLOT and ANALYSIS 
fig = figure;
semilogy(tspan, abs(x_analytical(tspan) - x_RK4), 'k.-')
hold on
semilogy(tspan, abs(x_analytical(tspan) - x_ABM3), '.-', 'Color', '#0011ee')
% semilogy(tspan, x_ABM3, '.-', 'Color', '#ee1100')
% plot(tspan, x, '.-', 'Color', '#11ee00')

hold off
xlabel('Time [s]')
ylabel('x value')
title('ODE solutions with variable step and method')
grid on
% Make it automatic with cells
legend('Analytical', ['RK2 Heun h = ', num2str(h(1))], ['RK4 h = ', num2str(h(1))]);


