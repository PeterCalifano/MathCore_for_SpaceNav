%% Modeling and Simulation of Aerospace Systems 2022/2023
% Assignment 1
% Author: Pietro Califano


%% Exercise 1
close all; 
clear; clc

DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

f = @(x) sin(x) + x - 1;

% Verify the function plot to check useful interval
a = 0; b = 1; % Range chosen by looking at the plot

figure(1);
x_plot  = linspace(-3, 3, 200);
y_plot = sin(x_plot) + x_plot - 1;

% Function plot
plot(x_plot, y_plot, '-', 'Color', '#0033dd', 'LineWidth', 1.17)
hold on
% Interval bounds plots
yline(0, 'LineWidth', 1.05)
xline(0, 'LineWidth', 1.05)
plot([a, a], [0 f(a)], 'k.-', 'MarkerSize', 15);
plot([b, b], [0 f(b)], 'k.-', 'MarkerSize', 15);

% Point labels
Fontsize = 15;
text(a-0.1, 0.09, '$\mathbf{a}$', 'FontSize', Fontsize)
text(b+0.04, 0.09, '$\mathbf{b}$', 'FontSize', Fontsize)
text(a+0.08, f(a), '$\mathbf{f(a)}$', 'FontSize', Fontsize)
text(b-0.1*b, f(b)+ 0.4*f(b), '$\mathbf{f(b)}$', 'FontSize', Fontsize)

% Plot options
hold off
ax = gca;

ax.XAxisLocation = 'origin'; 
ax.YAxisLocation = 'origin';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
xlabel('$\bf{x}$');
ylabel('$\bf{f(x)}$');
grid on

axis auto


tol = 1e-9;

% method = 1; % 1: Bi-section, 2: Secants, 3: Regula Falsi

n_samples = 100;
timer = nan(n_samples, 3);
x_zero = nan(1, 3);

% Call FindZero function to solve the problem
for method = 1:3
    for id = 1:n_samples
        tic;
        [x_zero(method), feval, ~] = FindZero(f, [a, b], method, tol);
        timer(id, method) = toc;
    end
end

figure(1)
hold on

plot(x_zero(method), 0, 'ko', 'Markersize', 3, 'LineWidth', 2)
text(x_zero(method)-0.10, -0.28, '$\mathbf{f(x)=0}$', 'FontSize', Fontsize, 'FontWeight','bold')

legend('$\mathbf{f(x) = sin(x)+x-1}$')

% Sample mean and std. dev. of computational time
Bisect_mean = mean(timer(:, 1));
Bisect_stddev = std(timer(:, 1));

Secants_mean = mean(timer(:, 2));
Secants_stddev = std(timer(:, 2));

RegFalsi_mean = mean(timer(:, 3));
RegFalsi_stddev = std(timer(:, 3));

check_estimate_convergence = 0;

if check_estimate_convergence == 1
    % Compute cumulative sample mean and standard dev.
    plot_samples = round(linspace(1, length(timer), 1000));

    for method = 1:3
        figure(method+1)
        hold on;
        cum_mean = nan();
        cum_stddev = nan();

        for cumid = 1:length(plot_samples)

            cum_mean(cumid) = mean(timer(1:plot_samples(cumid), method));
            cum_stddev(cumid) = std(timer(1:plot_samples(cumid), method));

        end

        subplot(2, 1, 1);
        plot(cum_mean, '.-');
        xlabel('Evaluation points')
        ylabel('Value')
        title('Sample $\hat{\mu}$ dynamics', 'Interpreter', 'latex')

        subplot(2, 1, 2);
        plot(cum_stddev, '.-');
        xlabel('Evaluation points')
        ylabel('Value')
        title('Sample $\hat{\sigma}$ dynamics', 'Interpreter', 'latex')
    end
end

%% Exercise 2
close all; 
clear; clc

DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

% Define function
syms x [2 1]
f = [x(1)^2 + x(2) - 5;
     x(2)^2 - x(1)];

cycle_iter = cell(2);
xsol = cell(2);
error = cell(2, 1);

% Analytical derivative
df = jacobian(f); % [2x2]

% Solution cycle
for guess = 1:2

    safe_var = 0;

    if guess == 1
        x_k = [1, 1]';
    elseif guess == 2
        x_k = [3, -2]';
    end

    err = [1, 1]';
    tol = 1e-8;

    while err(1) > tol || err(2) > tol
        if safe_var == 100
            break;
        end

        f_k = eval(subs(f, x, x_k));
        x_new = eval(x_k - subs(df, x, x_k)\f_k);

        f_xnew = eval(subs(f, x, x_new));

        x_k = x_new;
        err = abs(f_xnew);
        safe_var = safe_var + 1;

    end

    % Number of iterations of Jacobian
    cycle_iter{guess, 1} = safe_var;
    % Save solution from Jacobian
    xsol{guess, 1} = x_new;
    x_J = x_new;

    clear x_k x_new f_k 


    % FDM solution cycle

    safe_var = 0;

    if guess == 1
        x_k = [1, 1]';
    elseif guess == 2
        x_k = [3, -2]';
    end

    f_k = eval(subs(f, x, x_k));

    while abs(f_k(1)) > tol || abs(f_k(2)) > tol
        if safe_var == 100
            break;
        end

        df_k = nan(2);
        dx = zeros(length(f));

        for id = 1:length(f)
            % Define perturbation magnitude dx
            dx(id, id) = 0.001; %sqrt(eps)*max(1, abs(x_k(id)));
            f_xpert = eval(subs(f, x, x_k+dx(:, id)));
            df_k(:, id) = (f_xpert - f_k)./dx(id, id);
        end

        % Compute x value for next step
        x_new = x_k - df_k\f_k;

        % Store value at kth step
        f_kprev = f_k;

        % Update trial value
        f_k = eval(subs(f, x, x_new));
        x_k = x_new;
        safe_var = safe_var + 1;

    end
    % Number of iterations of FDM
    cycle_iter{guess, 2} = safe_var;
    % Save solution from FDM
    xsol{guess, 2} = x_new;
    % Compare it with Jacobian
    error{guess, 1} = abs(x_new - x_J);
end

% Relative error computation for second zero (guess = 2, row of the cell)
RE{1} = abs(xsol{2, 2} - xsol{2, 1})./abs(xsol{2, 1});


% PLOT SECTION

[X, Y] = meshgrid(-5:0.01:5, -5:0.01:5);

figure(1)

eval1 = X.^2 + Y - 5;
eval2 = Y.^2 - X;

% Fcn contours
contour(X, Y, eval1, [0, 0], 'Color', '#dd4411', 'LineWidth', 1.15);
hold on;
contour(X, Y, eval2, [0, 0], 'Color', '#44dd22', 'LineWidth', 1.15);

xlabel('$x_1$');
ylabel('$x_2$');
grid on

% Fcn zeros
plot(xsol{1, 1}(1), xsol{1, 1}(2), '.', 'Color', '#0022ee', 'MarkerSize', 15);
plot(xsol{2, 1}(1), xsol{2, 1}(2), '.', 'Color', '#0022ee', 'MarkerSize', 15);

% hold off
% ax.XAxisLocation = 'origin'; 
% ax.YAxisLocation = 'origin';
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XTick = -5:1:5;
ax.YTick = -5:1:5;
text(2.1, 1.15, '$\mathbf{z}_1 = [1.903, 1.379]$');
text(2.6, -1.4, '$\mathbf{z}_2 = [2.570, -1.603]$');

plot(1, 1, 'Marker', 'o', 'MarkerSize', 4, 'Color', 'k')
text(1.05, 0.85, '$\mathbf{z^1_{guess}}$')

plot(3, -2, 'Marker', 'o', 'MarkerSize', 4, 'Color', 'k')
text(3.05, -2.15, '$\mathbf{z^2_{guess}}$')
hold off
legend('Zeros of $f_{1} = {x_1}^2 + x_2 - 5$', 'Zeros of $f_{2} = {x_2}^2 - x_1$')

axis auto

%% Exercise 3
close all; 
clear; clc

rng default; % for repeatibility

DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)



% Options for statistical analysis
n_points = 200;
n_samples = 10000;

% Create non-uniformly spaced array
LowerB = 0.001;
UpperB = 0.4;
alfa = 0.999; % decides where and how much the point density is increased
csi = linspace(LowerB^(1/alfa), UpperB^(1/alfa), n_points); % uniformly spaced mesh to be passed to the mapping function
h_samples = csi.^alfa;

% Flag to switch last plot on or off
correlation_error_time = 1;

% Equation: xdot = -(t^2 - 1)*x/(t^2 + 1)
x0 = 1; % at time t0 = 0
x_analytical = @(t) exp(2*atan(t) - t); 

% Time span bounds
t0 = 0;
tf = 2;

% RHS of the system
f = @(x, t) -(t.^2 - 1).*x./(t.^2 + 1);
% Step size
h = [0.5, 0.2, 0.05, 0.01];

initialColorOrder = get(gca,'ColorOrder');
Colors = initialColorOrder(1:7, :);

% Variable preallocation
evalcounter_RK2 = nan(1, length(h));
evalcounter_RK4 = nan(1, length(h));
exetime = nan(n_samples, 4, 2);
x_ref = cell(4, 1);
err_RK2 = cell(4, 1);
err_RK4 = cell(4, 1);
time_vec = cell(4, 1);

figure(1);
% Evaluation with the assigned step-sizes
for idh = 1:4
    % Select h
    tspan = t0:h(idh):tf;
    x_ref{idh} = x_analytical(tspan);

    %         if tspan(end) < tf
    %             tspan = [tspan, tf];
    %         end

    for id_sample = 1:n_samples
       
        time_vec{idh} = tspan;

        % RK2 Heun
        exetime_clockRK2 = tic;
        [x_RK2, evalcounter_RK2(idh)] = RKN_integrator(f, tspan, x0, 2);
        exetime(id_sample, idh, 1) = toc(exetime_clockRK2);

        % Error evaluation for RK2 wrt analytical solution
        err_RK2{idh} = abs(x_RK2 - x_ref{idh}');

        % RK4
        exetime_clockRK4 = tic;
        [x_RK4, evalcounter_RK4(idh)] = RKN_integrator(f, tspan, x0, 4);
        exetime(id_sample, idh, 2) = toc(exetime_clockRK4);
        % Error evaluation for RK4 wrt analytical solution
        err_RK4{idh} = abs(x_RK4 - x_ref{idh}');
    end

    % PLOT and ANALYSIS
    subplot(2, 2, idh)
    analyticaltspan = t0:0.01:tf;
    plot(analyticaltspan, x_analytical(analyticaltspan), 'k-', 'LineWidth', 1.05)
    hold on
    plot(tspan, x_RK2, '.-', 'Color', Colors(1, :), 'LineWidth', 1.02)
    plot(tspan, x_RK4, '.-', 'Color', Colors(2, :), 'LineWidth', 1.02)
    % plot(tspan, x, '.-', 'Color', '#11ee00')

    hold off
    xlabel('$Time$')
    ylabel('$\mathbf{x(t)}$')
    grid minor

    title("Solutions with step h = " + num2str(h(idh)), 'FontWeight', 'bold', 'FontSize', 16)
    axis equal

    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    legend('Analytical', 'RK2 Heun', 'RK4 h');
end

% Pre-allocation of stats indices
mean_errRK2 = nan(4, 1);
mean_errRK4 = nan(4, 1);
mean_timeRK2 = nan(4, 1);
std_timeRK2 = nan(4, 1);
mean_timeRK4 = nan(4, 1);
std_timeRK4 = nan(4, 1);

% Statistics of errors
for idh = 1:4
    mean_errRK2(idh) = mean(err_RK2{idh});
    mean_errRK4(idh) = mean(err_RK4{idh});

    mean_timeRK2(idh) = mean(exetime(:, idh, 1));
    std_timeRK2(idh) = std(exetime(:, idh, 1));
    mean_timeRK4(idh) = mean(exetime(:, idh, 2));
    std_timeRK4(idh) = std(exetime(:, idh, 2));

end

% legend('Analytical', ['RK2 Heun h =', num2str(h(idh))], ['RK4 h = ', num2str(h(idh))]);

% First plot: Integration error vs time t
RK2_err_fig = figure;
RK4_err_fig = figure;
markerlist = ['v', 'd', 's', '^'];

for idh = 1:length(h)

    figure(RK2_err_fig);
    semilogy(time_vec{idh}, err_RK2{idh}, '-', 'Marker', markerlist(idh));
    hold on

    figure(RK4_err_fig);
    semilogy(time_vec{idh}, err_RK4{idh}, '-', 'Marker', markerlist(idh));
    hold on
    
end

legend_cell= {"h = " + num2str(h(1)), "h = " + num2str(h(2)), "h = " + num2str(h(3)), "h = " + num2str(h(4))};

figure(RK2_err_fig);
xlabel('Time')
ylabel('$|\bf{x}_{RK2} - \bf{x}_{analytical}|$')
grid on

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
legend(legend_cell, 'Location', 'southeast')

figure(RK4_err_fig);
xlabel('Time')
ylabel('$|\bf{x}_{RK4} - \bf{x}_{analytical}|$')
grid on

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
legend(legend_cell, 'Location', 'southeast')

% Computational time vs Integration error plot
Meantime = nan(length(h), 2);
StdDevtime = nan(length(h), 2);

for idh = 1:4

    Meantime(idh, 1) = mean(exetime(:, idh, 1));
    Meantime(idh, 2) = mean(exetime(:, idh, 2));

    StdDevtime(idh, 1) = std(exetime(:, idh, 1));
    StdDevtime(idh, 2) = std(exetime(:, idh, 2));

end

pause(0.001)
if correlation_error_time == 1
warning(['Analysis of integration error vs time started.' ...
    ' This may take a while to run!']);

% Statistical analysis cycle
% Vectors of [n_points x 2]: 1st column: Mean int.error. 2nd column:
point_samples_RK2 = nan(n_points, 2);
point_samples_RK4 = nan(n_points, 2);

tic
for idh = 1:length(h_samples)

    % Select h and build tspan
    tspan = t0:h_samples(idh):tf;
    % Compute analytical solution for reference
    x_ref = x_analytical(tspan);

    exetime_statRK2 = nan(n_samples, 1);
    exetime_statRK4 = nan(n_samples, 1);

    % Execute multiple times to estimate computational time
    for id_sample = 1:n_samples

        % RK2 Heun
        exetime_clockRK2 = tic;
        [x_RK2, ~] = RKN_integrator(f, tspan, x0, 2);
        exetime_statRK2(id_sample, 1) = toc(exetime_clockRK2);

        % RK4
        exetime_clockRK4 = tic;
        [x_RK4, ~] = RKN_integrator(f, tspan, x0, 4);
        exetime_statRK4(id_sample, 1) = toc(exetime_clockRK4);

    end

    % Average computational time RK2
    point_samples_RK2(idh, 1) = mean(exetime_statRK2);
    % Error evaluation for RK2 wrt analytical solution
    point_samples_RK2(idh, 2) = mean(abs(x_RK2 - x_ref'));

    % Average computational time RK4
    point_samples_RK4(idh, 1) = mean(exetime_statRK4);
    % Error evaluation for RK4 wrt analytical solution
    point_samples_RK4(idh, 2) = mean(abs(x_RK4 - x_ref'));
end


stat_analysis = figure;
loglog(point_samples_RK2(:, 2), point_samples_RK2(:, 1), '.','Color', '#bb5511', 'MarkerSize', 10)
hold on
loglog(point_samples_RK4(:, 2), point_samples_RK4(:, 1), '.','Color', '#2244bb', 'MarkerSize', 10)

grid on
xlabel('Mean integration error ($|x_{num}(t) - x_{ref}(t)|$)')
ylabel('Computational time [s]')
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
legend('RK2', 'RK4');

end
%% Exercise 4
close all;
clear; clc

DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)
color_palette = {'#ee2200', '#00ee22', '#ffdc20', '#0022ee'};

% Operators definitions
A = @(a) [0 1;
    -1 2*cos(a)];

F_RK2Heun = @(h, a) eye(2) + h*A(a) + 0.5*(h^2 .* A(a)^2);
F_RK4 = @(h, a) eye(2) + h.*(A(a) + 0.5.*h.*(A(a)^2) + (1/6)*h.^2*(A(a)^3) + (1/24) *h^3*A(a)^4);

scheme_selector = {F_RK2Heun, F_RK4};

% Options
n_points = 180*3;
alpha_range = deg2rad(linspace(180, 0, n_points));
x0 = [4, 4];
legend_cell = {'RK2 Heun', 'RK4'};

% Stability domain boundary plot
for idscheme = 1:numel(scheme_selector)
    guess = x0(idscheme);
    % Call to PlotRos function (solves the problem to find the boundary and plots it)
    [RoS_bounds, h_max] = PlotRoS(scheme_selector{idscheme}, alpha_range, guess, 0, 1);

    if idscheme == 1

        initialColorOrder = get(gca,'ColorOrder');
        Colors = initialColorOrder(1:2, :);

    end

    patch('Xdata', real(RoS_bounds), 'Ydata', imag(RoS_bounds), 'FaceColor', Colors(idscheme, :), 'FaceAlpha', 0.2)

    if idscheme == 1
        RoS_bounds_RK2 = [RoS_bounds, flip(conj(RoS_bounds))];
    elseif idscheme == 2
        RoS_bounds_RK4 = [RoS_bounds, flip(conj(RoS_bounds))];
    end
end

% hold off;

yline(0, 'LineWidth', 1.05)
xline(0, 'LineWidth', 1.05)
plot(-2.78528, 0, 'k.', 'MarkerSize', 14)
text(-3.65, 0.1, '$\mathbf{h\lambda_{RK4}} = -2.78$');
plot(-2, 0, 'k.', 'MarkerSize', 14)
text(-1.95, -0.1, '$\mathbf{h\lambda_{RK2}} = -2$');
hold off;

legend('RK2 RoS boundary', 'RK2 RoS', 'RK4 RoS boundary', 'RK4 RoS', 'Location', 'best')
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

% legend(legend_cell);

% Point 2)
alpha = pi;
x0 = 4;
F = F_RK2Heun;
fun = @(h, a) max(abs(eig(F(h, a)))) - 1;
opts = optimoptions('fsolve', 'FunctionTolerance', 1e-12, 'Display', 'off', 'UseParallel', false);

[h_RK2_pi, ~, exitflag_hRK2pi, ~] = fsolve(@(h) fun(h, alpha), x0, opts);

% Point 3)
% Time span bounds
t0 = 0;
tf = 2;

% RHS of the system
% f = @(x, t) -(t^2 - 1)*x./(t^2 + 1);

% Eigenvalue of the system (scalar)
lambda = @(t) -(t.^2 - 1)./(t.^2 + 1);
% Step size
h = [0.5, 0.2, 0.05, 0.01];

hlambda_cell = cell(4, 1);

figure(2);
% Outer For cycle: step size h
for idh = 1:length(h)
    % Define time span
    tspan = t0:h(idh):tf;

%     if tspan(end) < tf
%         tspan(end+1) = tf;
%     end

    % Inner For cycle: time t
    hlambda_intime = h(idh).*lambda(tspan);

    % Save hlambdas
    hlambda_cell{idh} = hlambda_intime;

    % Zoom to show the eigenvalues in last point
    hold on
    scatter(hlambda_cell{idh}, zeros(1, length(hlambda_cell{idh})), 'x', 'Color', color_palette{idh}, 'Linewidth', 1);

    legend_cell{idh} = "h = " + num2str(h(idh));
    clear hlambda_intime

end

% Plot options
legend_cell{idh+1} = "RK4 RoS boundary";
legend_cell{idh+2} = "RK2 RoS boundary";


plot(real(RoS_bounds_RK4), imag(RoS_bounds_RK4), 'k-', 'Linewidth', 1.1);
hold on
plot(real(RoS_bounds_RK2), imag(RoS_bounds_RK2), '-', 'Linewidth', 1.1, 'Color', '#aa0000');

patch('Xdata', real(RoS_bounds_RK4), 'Ydata', imag(RoS_bounds_RK4), 'FaceColor', '#000000', 'FaceAlpha', 0.1)

patch('Xdata', real(RoS_bounds_RK2), 'Ydata', imag(RoS_bounds_RK2), 'FaceColor', '#aa0000', 'FaceAlpha', 0.2)

xlabel('Re$\{h\lambda\}$')
ylabel('Im$\{h\lambda\}$')

axis equal
grid on

ax = gca;

ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
yline(0, 'LineWidth', 1.05)
xline(0, 'LineWidth', 1.05)

legend(legend_cell)
hold off;

color_palette = {'#ee2200', '#00ee22', '#ffdc20', '#0022ee'};


% figure(3)
% Second reference frame in box
ax2 = axes('Position', [0.4 0.25 0.5 0.5], 'Box','on');
hold on;

for idh = 1:4
    scatter(hlambda_cell{idh}, zeros(1, length(hlambda_cell{idh})), 'x',  'Color', color_palette{idh}, 'Linewidth', 1);
end

plot(real(RoS_bounds_RK4), imag(RoS_bounds_RK4), 'k-', 'Linewidth', 1.1);
hold on
plot(real(RoS_bounds_RK2), imag(RoS_bounds_RK2), '-', 'Linewidth', 1.1, 'Color', '#aa0000');
patch('Xdata', real(RoS_bounds_RK4), 'Ydata', imag(RoS_bounds_RK4), 'FaceColor', 'k', 'FaceAlpha', 0.1)
patch('Xdata', real(RoS_bounds_RK2), 'Ydata', imag(RoS_bounds_RK2), 'FaceColor', '#aa0000', 'FaceAlpha', 0.2)
yline(0, 'LineWidth', 1.05)
xline(0, 'LineWidth', 1.05)
% legend(legend_cell)
xlabel('Re$\{h\lambda\}$')
ylabel('Im$\{h\lambda\}$')
axis equal

ax2.XMinorTick = 'on';
ax2.YMinorTick = 'on';
grid on

xlim([-0.06 0.02]);
ylim([-0.1 0.1]);

% legend(legend_cell)

%% Exercise 5
close all; 
clear; clc

DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

% Parameters
x0 = [1; 1]; % Initial Condition
alpha_range = deg2rad(linspace(180, 0, 100));
t0 = 0;
tf = 1;
color_palette = {'#ee2200', '#00ee22', '#eece00', '#0022ee'};

% Functions definition
A = @(a) [0 1;
    -1 2*cos(a)];

F_FE = @(h, a) eye(2) + h.*A(a);
F_RK2Heun = @(h, a) eye(2) + h*A(a) + 0.5*(h^2 .* A(a)^2);
F_RK4 = @(h, a) eye(2) + h.*(A(a) + 0.5.*h.*(A(a)^2) + (1/6)*h.^2*(A(a)^3) + (1/24) *h^3*A(a)^4); 

x_analytical = @(a, time) expm(A(a)*time)*x0;

% Error at final time
tol = [1e-3, 1e-4, 1e-5, 1e-6];

% Very sensitive to initial guess. Iterate a bit also looking at the
% function plot 

% h_guess = [1e-3, 0.001, 0.0005, 0.00005;
%     0.08, 0.04, 0.01, 0.005;
%     0.055, 0.04, 0.004, 0.1];

h_guess = [1e-4, 1e-4, 1e-5, 1e-6;
    0.048, 0.015, 0.019, 0.008;
    0.551, 0.283, 0.15, 0.06];


% Time vector for the numerical solution
% opts = optimoptions('fsolve', 'FunctionTolerance', 1e-12, 'Display', 'off', 'UseParallel', false);
opts = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'StepTolerance', 1e-12,...
    'FunctionTolerance', 1e-14, 'Display', 'off', 'UseParallel', false);

Operator_list = {F_FE, F_RK2Heun, F_RK4};
Title_list = {'RK1', 'RK2 Heun', 'RK4'};

h_tol_data = nan(length(alpha_range), length(h_guess), length(Operator_list));
Tol_bounds = nan(2*length(alpha_range), 4, 3);

% Solve problem for each scheme and for each tolerance
for scheme_id = 1:length(Operator_list)
    op = Operator_list{scheme_id};

    for idh = 1:length(h_guess(scheme_id, :))
        
        for alpha_id = 1:length(alpha_range) % Solve for each alpha

            F = @(h) op(h, alpha_range(alpha_id)); 

            % Reference solution at final time
            x_ref = x_analytical(alpha_range(alpha_id), tf);

            % Zero-finding problem 
            zerofun = @(h) norm(x_ref - RK_evalfinalstate(F, x0, h, t0, tf), "inf") - tol(idh);

%           [h_tol, feval, exitflag, output] = fsolve(@(h) zerofun(h), h_guess(scheme_id, idh), opts);
            [h_tol, feval, exitflag, output] = lsqnonlin(@(h) zerofun(h), h_guess(scheme_id, idh), 0, 10, opts);

            h_tol_data(alpha_id, idh, scheme_id) = h_tol;
            %         feval_data(alpha_id, idh) = feval;
            %         exitflag_data(alpha_id, idh) = exitflag;
            %         output_data{alpha_id, idh} = output;

        end
        disp(['Cycle with to = ' num2str(tol(idh)), ' of scheme ', Title_list{scheme_id} ,' completed']);
    end
   
    % Assign h
    lambda = h_tol_data(:, :, scheme_id) .*repmat(((cos(alpha_range) + 1j* sin(alpha_range)))', 1, 4);
    % Build curve of the solution
    Tol_bounds_half = lambda;
    Tol_bounds(:, :, scheme_id) = [Tol_bounds_half; flip(conj(Tol_bounds_half))];


    for idh = 1:length(h_guess(scheme_id, :))
        figure(scheme_id)

        if ~exist('ax1', 'var')
            ax1 = gca;
        end

        hold(ax1);
        plot(ax1, real(Tol_bounds(:, idh, scheme_id)), imag(Tol_bounds(:, idh, scheme_id)), '-', 'LineWidth', 1.02, 'Color', color_palette{idh})
        hold on

        % Plot options
        xlabel('$Re\{h\lambda\}$')
        ylabel('$Im\{h\lambda\}$')
        axis equal
        grid on

        ax1.XAxisLocation = 'bottom';
        ax1.YAxisLocation = 'left';
        ax1.XMinorTick = 'on';
        ax1.YMinorTick = 'on';
        hold off;

        if scheme_id == 1 && (idh == 3 || idh == 4)
            if ~exist('ax2', 'var')
                ax2 = axes('Position', [0.68 0.38 0.26 0.26]);
            end
            hold(ax2)
            plot(ax2, real(Tol_bounds(:, idh, scheme_id)), imag(Tol_bounds(:, idh, scheme_id)), '-', 'LineWidth', 1.02, 'Color', color_palette{idh})
            hold on
            axis equal
            grid on
            ax2.XMinorTick = 'on';
            ax2.YMinorTick = 'on';

            hold off;
        end
    end
    legend(ax1, 'tol = 1e-3', 'tol = 1e-4', 'tol = 1e-5', 'tol = 1e-6');
    clear ax1 ax2
    disp(' ')

end

alpha_id = find(alpha_range == pi); % Extract id corresponding to pi in alpha_range
req_evals = [1, 2, 4]';

func_evals = nan(4 ,3);

figure(scheme_id+1);
for scheme_id = 1:3

    for idh = 1:numel(tol)

        nsteps = floor((tf-t0)/h_tol_data(alpha_id, idh, scheme_id));

        if tf - h_tol_data(alpha_id, idh, scheme_id)*nsteps > 0
            nsteps = nsteps + 1;
        end

        func_evals(idh, scheme_id) = nsteps * req_evals(scheme_id);
    end        % LOGLOG PLOT
    loglog(tol', func_evals(:, scheme_id), '-', 'Color', color_palette{scheme_id}, 'LineWidth', 1.1);
    hold on
end

hold off
legend(Title_list, 'Location', 'Best')

grid on
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
xlabel('Tolerance')
ylabel("Function evaluations")
% title('Function evals vs Tolerance')


%% Exercise 6
close all; 
clear; clc

DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)


A = @(a) [0 1;
    -1 2*cos(a)];

F_RK2Heun = @(h, a) eye(2) + h.*A(a) + 0.5*(h.^2 .* A(a)^2);

opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

x0 = 2.5;
alfa_range = deg2rad(0:0.1:180);

for theta = [0.1, 0.3, 0.4, 0.7, 0.9]

    F_BI1_theta = @(h, a) (eye(2) - (1-theta)*h*A(a))\(eye(2) + theta*h*A(a)); % FE + BE, validated
    % F_BI2_theta = @(h, a) (eye(2) + (theta-1)*h*A(a) + 0.5*(theta-1)^2*h^2*A(a)^2)\(eye(2) + theta*h*A(a) + 0.5*theta^2*h^2*A(a)^2) ;% BRK2: Heun forward and Heun backward
    F_BI2_theta = @(h, a) F_RK2Heun((theta-1)*h, a)\F_RK2Heun(theta*h, a);

    F_list = {F_BI1_theta, F_BI2_theta};
    F = F_list{2};

    fun = @(h, a) max(abs(eig(F(h, a)))) - 1;
    RoS_bound_temp = nan(1, length(alfa_range));

    % Find RoS boundary
    for id = 1:length(alfa_range)

        [h_max(id), feval, exitflag, output] = fsolve(@(h) fun(h, alfa_range(id)), x0, opts);
        RoS_bound_temp(id) = h_max(id)*(cos(alfa_range(id)) + 1j*sin(alfa_range(id)));  

    end

    RoS_bound = [RoS_bound_temp, conj(RoS_bound_temp)];

    figure(1)
    hold on
    scatter(real(RoS_bound), imag(RoS_bound), '.')
    clear RoS_bound

end

% Plot options
figure(1)
hold off

grid on

xlabel('Re$\{h\lambda\}$')
ylabel('Im$\{h\lambda\}$')

ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

yline(0, 'Color', 'k', 'LineWidth', 1.05)
xline(0, 'Color', 'k', 'LineWidth', 1.05)

initialColorOrder = get(gca,'ColorOrder');

Colors = initialColorOrder(1:5, :);

text(9.2, -4, '$\mathbf{BI2_{0.4}}$', 'Color', Colors(3, :), 'FontSize', 14, 'FontWeight', 'bold');
text(2.4, 0.4, '$\mathbf{BI2_{0.1}}$', 'Color', Colors(1, :), 'FontSize', 14, 'FontWeight', 'bold');
text(5.03, 1, '$\mathbf{BI2_{0.3}}$', 'Color', Colors(2, :), 'FontSize', 14, 'FontWeight', 'bold');
text(-6, -3, '$\mathbf{BI2_{0.7}}$', 'Color', Colors(4, :), 'FontSize', 14, 'FontWeight', 'bold');
text(-2.5, 0.4, '$\mathbf{BI2_{0.9}}$', 'Color', Colors(5, :), 'FontSize', 14, 'FontWeight', 'bold');

% legend_cell = {'$BI2_{0.1}$', '$BI2_{0.3}$', '$BI2_{0.4}$', '$BI2_{0.7}$', '$BI2_{0.9}$'};
% legend(legend_cell)

axis equal
% axis([-7 12 -6 6])
xlim([-6 12])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 7
close all; 
clear; clc

DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

initialColorOrder = get(gca,'ColorOrder');
Colors = initialColorOrder(1:5, :);

B = [-180.5, 219.5;
    179.5, -220.5];

% A = @(a) [0 1;
%     -1 cos(a)];

RHS_IEX4 = @(~, x) B*x;
RHS_RK4 = @(x, ~) B*x;

% Analytical solution
t0 = 0;
tf = 5;

h = 0.1;

x0 = [1; 1];
sol = @(t) expm(B.*t)*x0;

tspan = t0:h:tf;
tspan_plot = t0:0.01:tf;

if tspan(end) < tf
    tspan(end+1) = tf;
end

x_ref = nan(length(tspan), length(x0));
x_ref_plot = nan(length(tspan_plot), length(x0));

for idt = 1:length(tspan)
    x_ref(idt, :) = sol(tspan(idt));
end

for idt = 1:length(tspan_plot)
    x_ref_plot(idt, :) = sol(tspan_plot(idt));
end

% Solution with RK4
[x_RK4, ~] = RKN_integrator(RHS_RK4, tspan, x0, 4);

% Solution with IEX4
[x_IEX4, evalcounter_IEX4] = IEX4_integrator(RHS_IEX4, tspan, x0);

% IEX4 Stability domain analysis
A = @(a) [0 1;
        -1 2*cos(a)];
F_IEX4 = @(h, a) -1/6 * (eye(2) - h*A(a))^-1 + 4 * (eye(2) - 0.5 * h * A(a))^-2 - 27/2 * (eye(2) - 1/3 * h*A(a))^-3 + 32/3 * (eye(2) - 1/4 * h * A(a))^-4;
alpha_range = deg2rad(180:-1:0);
IEX4_guess = 5.05;

% Plot stability region
[RoS_bounds_IEX4, h_max_IEX4] = PlotRoS(F_IEX4, alpha_range, IEX4_guess, 0, 1);

patch('XData', real(RoS_bounds_IEX4), 'YData', imag(RoS_bounds_IEX4), 'FaceColor', Colors(1, :), 'FaceAlpha', 0.3, 'EdgeColor',  Colors(1, :));
% hold off;

F_RK4 = @(h, a) eye(2) + h.*(A(a) + 0.5.*h.*(A(a)^2) + (1/6)*h.^2*(A(a)^3) + (1/24) *h^3*A(a)^4);
RK4_guess = 4;

[RoS_bounds_RK4, h_max_RK4] = PlotRoS(F_RK4, alpha_range, RK4_guess, 0, 1);
patch('XData', real(RoS_bounds_RK4), 'YData', imag(RoS_bounds_RK4), 'FaceColor', Colors(2, :), 'FaceAlpha', 0.3, 'EdgeColor',  Colors(2, :));

% System eigenvalues (not time-varying)
eig_B = eig(B);
plot(h*eig_B(1), 0, 'x', 'Color', '#ee2222', 'MarkerSize', 10, 'Linewidth', 1.1);
plot(h*eig_B(2), 0, 'x', 'Color', '#ee2222', 'MarkerSize', 10, 'Linewidth', 1.1);

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
yline(0, 'LineWidth', 1.05)
xline(0, 'LineWidth', 1.05)

legend('IEX4: RoS boundary', 'IEX4: \textbf{Unstable} region',...
    "RK4: RoS boundary", "RK4: \textbf{Stable} region", "$h\lambda_1 = (-0.1,0)$", "$h\lambda_2 = (-40,0)$")

ax2 = axes('Position', [0.3 0.25 0.2 0.2], 'Box', 'on');
plot(real(RoS_bounds_RK4), imag(RoS_bounds_RK4), 'Color', Colors(2, :), 'LineWidth', 1.05)
patch('XData', real(RoS_bounds_IEX4), 'YData', imag(RoS_bounds_IEX4), 'FaceColor', Colors(1, :), 'FaceAlpha', 0.3, 'EdgeColor',  Colors(1, :));
hold on
plot(real(RoS_bounds_IEX4), imag(RoS_bounds_IEX4), 'Color', Colors(1, :), 'LineWidth', 1.05)
patch('XData', real(RoS_bounds_RK4), 'YData', imag(RoS_bounds_RK4), 'FaceColor', Colors(2, :), 'FaceAlpha', 0.3, 'EdgeColor',  Colors(2, :));

plot(h*eig_B(1), 0, 'x', 'Color', '#ee2222', 'MarkerSize', 10, 'Linewidth', 1.1);
plot(h*eig_B(2), 0, 'x', 'Color', '#ee2222', 'MarkerSize', 10, 'Linewidth', 1.1);

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
yline(0, 'LineWidth', 1.05)
xline(0, 'LineWidth', 1.05)
axis equal
grid on
xlim([-0.4 0.4])
ylim([-0.5 0.5])

% Error computation
err_RK4 = abs(x_ref - x_RK4);
err_IEX4 = abs(x_ref - x_IEX4);

% Temporary plots
% Solutions
figure(2)
subplot(2, 1, 1)
plot(tspan_plot, x_ref_plot(:, 1), '-', 'Color', '#1144bb', 'Linewidth', 1.1)
hold on
plot(tspan, x_IEX4(:, 1), '.', 'Color', '#ee4411', 'Markersize', 12)

xlabel('Time')
ylabel('$\mathbf{x_1(t)}$')
grid on
axis equal
legend('Analytical', 'IEX4')
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

subplot(2, 1, 2)
plot(tspan_plot, x_ref_plot(:, 2), '-', 'Color', '#1144bb' ,'Linewidth', 1.1)
hold on
plot(tspan, x_IEX4(:, 2), '.', 'Color', '#ee4411', 'Markersize', 12)

xlabel('Time')
ylabel('$\mathbf{x_2(t)}$')
grid on
axis equal
legend('Analytical', 'IEX4')
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

% Error analysis vs time
figure(3)
% scatter(tspan, err_RK4, 'Color', '#0011ee')
semilogy(tspan(2:end), err_IEX4(2:end, 1), '.-', 'Color', '#cc3300', 'LineWidth', 1.05)
hold on
semilogy(tspan(2:end), err_IEX4(2:end, 2), '.-', 'Color', '#0033cc', 'LineWidth', 1.05)


legend('$\mathbf{x_1}$ \textbf{err - IEX4}', '$\mathbf{x_2}$ \textbf{err - IEX4}')
xlabel('Time')
ylabel('$|\mathbf{x}_{sol} - \mathbf{x}_{IEX4}|$')
grid on
axis auto


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 8
close all; 
clear; clc

% NOTE: Plots of trajectory and integration errors have selectors before
% for cycle. ID legend: 1) RK4, 2) AB3, 3) AM3, 4) ABM3, 5) BDF3

DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

initialColorOrder = get(gca,'ColorOrder');
Colors = initialColorOrder(1:6, :);
markerlist = ['*', 'v', 'd', 's', '^'];
% Note: x = [x; y];
% Options
t0 = 0;
tf = 3;
h = 0.1; % Note h = 0.01 makes all the integrators converging to the same results

tspan = t0:h:tf;

if tspan(end) < tf
    tspan(end+1) = tf;
end

x0 = [1; 1];
RHS = @(x, t) [-5/2 * (1+8*sin(t)) * x(1);
    (1 - x(1))*x(2) + x(1)];

% Reference solution with ode45 and ode113
opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
[time_ode45, state_ode45] = ode45(@(t, x) RHS(x, t), t0:0.001:tf, x0, opts);

[time_ode113, state_ode113] = ode113(@(t, x) RHS(x, t), t0:0.001:tf, x0, opts);
state_ode113_resampled(:, 1) = interp1(time_ode113, state_ode113(:, 1), tspan, 'spline');
state_ode113_resampled(:, 2) = interp1(time_ode113, state_ode113(:, 2), tspan, 'spline');

% RK4
tic
[x_RK4, evalcounter_RK4] = RKN_integrator(RHS, tspan, x0, 4);
int_err{1, 1} = abs(state_ode113_resampled(:, 1) - x_RK4(:, 1));
int_err{1, 2} = abs(state_ode113_resampled(:, 2) - x_RK4(:, 2));
toc

% AB3
tic
[~, x_AB3, evalcounter_AB3] = Multistep_integrator(RHS, tspan, x0, 'AB', 3);
int_err{2, 1} = abs(state_ode113_resampled(:, 1) - x_AB3(:, 1));
int_err{2, 2} = abs(state_ode113_resampled(:, 2) - x_AB3(:, 2));
toc

% AM3
tic
[~, x_AM3, evalcounter_AM3] = Multistep_integrator(RHS, tspan, x0, 'AM', 3);
int_err{3, 1} = abs(state_ode113_resampled(:, 1) - x_AM3(:, 1));
int_err{3, 2} = abs(state_ode113_resampled(:, 2) - x_AM3(:, 2));
toc

% ABM3 
tic
[~, x_ABM3, evalcounter_ABM3] = Multistep_integrator(RHS, tspan, x0, 'ABM', 3);
int_err{4, 1} = abs(state_ode113_resampled(:, 1) - x_ABM3(:, 1));
int_err{4, 2} = abs(state_ode113_resampled(:, 2) - x_ABM3(:, 2));
toc

% BDF3
tic
[~, x_BDF3, evalcounter_BDF3] = Multistep_integrator(RHS, tspan, x0, 'BDF', 3);
int_err{5, 1} = abs(state_ode113_resampled(:, 1) - x_BDF3(:, 1));
int_err{5, 2} = abs(state_ode113_resampled(:, 2) - x_BDF3(:, 2));
toc

err_final = cell(5, 2);

for i = 1:5
    err_final{i, 1} = int_err{i, 1}(end);
    err_final{i, 2} = int_err{i, 2}(end);
end

figure(1)

AB3flag = 1; % Plot AB3 or the other schemes solutions

for i = 1:2

    subplot(2, 1, i);
    plot(tspan, state_ode113_resampled(:, i), 'k-', 'Linewidth', 1.05)
    hold on


    if AB3flag == 1
        plot(tspan, x_AB3(:, i), '-', 'Color',  Colors(2, :), 'Linewidth', 1.05)
    else
        %     plot(tspan, x_RK4(:, i), '-', 'Color',  Colors(4, :), 'Linewidth', 1.05)
        plot(tspan, x_AM3(:, i), '-', 'Color', Colors(1, :), 'Linewidth', 1.05)
        plot(tspan, x_ABM3(:, i), '-', 'Color',  Colors(2, :), 'Linewidth', 1.05)
        plot(tspan, x_BDF3(:, i), '-', 'Color',  Colors(3, :), 'Linewidth', 1.05)
    end

    grid on
    xlabel('Time')

    if i == 1
        ylabel('$\mathbf{x(t)}$', 'FontWeight', 'bold')
    elseif i == 2
        ylabel('$\mathbf{y(t)}$', 'FontWeight', 'bold')
    end

    if AB3flag == 1
        legend('ode113 solution', 'AB3', 'Location', 'best');
    else
        legend('ode113 solution','AM3', 'ABM3', 'BDF3', 'Location', 'best');
    end

    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    axis auto
end

% Integration error plot
% ID legend: 1) RK4, 2) AB3, 3) AM3, 4) ABM3, 5) BDF3
legend_cell_fixed = {'RK4', 'AB3', 'AM3', 'ABM3', 'BDF3'};

schemes = [3, 4, 5];
legendcell = cell(1, length(schemes));

figure(2)
for i = 1:2
    subplot(2, 1, i);
    hold on

    id = 1;
    for idscheme = schemes

        if i == 1
            semilogy(tspan, int_err{idscheme, 1}, 'Color', Colors(idscheme, :), 'LineWidth', 1.02, 'Marker', markerlist(idscheme));
        elseif i == 2
            semilogy(tspan, int_err{idscheme, 2}, 'Color', Colors(idscheme, :), 'LineWidth', 1.02, 'Marker', markerlist(idscheme));
        end

        legendcell{id} = legend_cell_fixed{idscheme};
        id = id + 1;
    end

    xlabel('Time', 'FontWeight', 'bold')

    if i == 1
        ylabel('$\mathbf{|x_{ode113}(t) - x(t)|}$')
    elseif i == 2
        ylabel('$\mathbf{|y_{ode113}(t) - y(t)|}$')
    end

    legend(legendcell, 'Location', 'best');
    grid on

    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    axis auto
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
% Not exercise specific

function [RoS_bounds, h_max] = PlotRoS(F, alpha_range, x0, half, plot_flag)
%% PROTOTYPE
% [RoS_bounds, h_max] = PlotRoS(F, alpha_range, x0, half, plot_flag)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Solve the problem "Find h >= 0 s.t. max(abs(eig(F(h, a)))) - 1 = 0" to
% generate the Region of stability boundary of Single step methods, specified by F operator,
% for a given range of alpha
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% F: [function handle] Operator as function of (h, alpha) (in this order)
% alpha_range: [1xN] range of alpha for which the problem must be solved
% x0: [1] initial guess for h stepsize
% half: [flag] 0 or 1 to obtain the entire boundary or only half
% plot_flag: [flag] 0 or 1 to switch the plot on or off
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% RoS_bounds: [1x2*N] complex numbers representing the RoS boundary
% h_max_ [1xN] values of h corresponding to the RoS boundary
% -------------------------------------------------------------------------------------------------------------

fsolve_tol = 1e-8;

if ~exist("plot_flag", 'var')
    plot_flag = 0;
end

%% Function code
% Plot settings
DefaultFontSize = 13;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');  
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

if iscell(F)
    scheme_n = length(F);
    warning('Plot management with cell input not completed');
else
    scheme_n = 1;
end

% Function options
opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', fsolve_tol);

% Functional for the zero-finding procedure
fun = @(h, a) max(abs(eig(F(h, a)))) - 1;

% Static memory allocation
h_max = nan(length(alpha_range), 1);
exitflag = nan(length(alpha_range), 1);

 for id = 1:length(alpha_range)
     [h_max(id), ~, exitflag(id), ~] = fsolve(@(h) fun(h, alpha_range(id)), x0, opts);
 end

 h_max(exitflag ~= 1) = [];
 alpha_range(exitflag ~= 1) = [];

 if scheme_n == 1
     RoS_bounds = h_max'.*(cos(alpha_range) + 1j*sin(alpha_range));
     if half ~= 1
         RoS_bounds = [RoS_bounds, flip(conj(RoS_bounds(2:end)))];
     end
 else
%      RoS_bounds{scheme_n} = h_max'.*(cos(alpha_range) + 1j*sin(alpha_range));
%      if half ~= 1
%          RoS_bounds{scheme_n} = [RoS_bounds{scheme_n}, conj(RoS_bounds{scheme_n})];
%      end
 end

 figure(1)
 hold on

 % Plot unit circle
%  phase = linspace(0, 2*pi, 100);
%  r = 1;
%  unit_circle = r*exp(1i*phase);
%  plot(real(unit_circle), imag(unit_circle), '-', 'Color', '#bb1111', 'LineWidth', 1.1);

if scheme_n == 1
    if plot_flag == 1
        plot(real(RoS_bounds), imag(RoS_bounds), '-', 'LineWidth', 1.10)
    elseif plot_flag == 0
        scatter(real(RoS_bounds), imag(RoS_bounds), '.');
    end
else
%     for id_scheme = 1:length(RoS_bounds)
%         scatter(real(RoS_bounds{id_scheme}), imag(RoS_bounds{id_scheme}), '.')
%     end
end


 axis equal
 grid on

 % Plot options
 xlabel('Re$\{h\lambda\}$')
 ylabel('Im$\{h\lambda\}$')
 title('Region of numerical stability - $\{h\lambda\}$ plane')

 % Plot axes
%  yline(0, 'LineWidth', 1.1);
%  xline(0, 'LineWidth', 1.1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCISE 1
function [x_zero, feval, info] = FindZero(fun, interval, method, tol)
%% PROTOTYPE
% [x_zero, feval, info] = FindZero(fun, interval, method, tol)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Analogous of fzero implementing Bisection, Secants and Regula Falsi
% algorithms to find zeros of a scalar function 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% fun: [function handle] any scalar function
% interval: [1x2] search interval specified as [a, b]
% method: [char] 1,2,3 = Bisection, Secants, Regula Falsi
% tol: [scalar] acceptable tolerance on the solution
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% x_zero: [scalar] value of the zero 
% feval: [scalar] value of the function at the found zero
% info: [struct] with field "evalcounter" counting function evals
% -------------------------------------------------------------------------------------------------------------

if nargin < 4
    tol = 1e-8;
    disp(['Default tolerance: ', num2str(tol)])
end

% Initialization of counters and variables
evalcounter = 0;
MaxIter = 200;
safety_var = 0;

a = interval(1);
b = interval(2);

f_a = fun(a);
evalcounter = evalcounter + 1;
f_b = fun(b);
evalcounter = evalcounter + 1;


if b < a
    warning('Left value of the interval seems to be greater than the right value. Zero finding procedure may fail.')
end

switch method

    case 1 % Bi-section
        c = b;
        err = 1;
        while abs(err) > tol
            if safety_var == MaxIter
                warning('Solution not found in specified Max Iteration number')
                break;
            end

            if f_a * f_b > 0
                error('Zero not existent in the selected interval');
            end

            c_prev = c;
            % Compute middle point of the interval
            c = (a+b)/2;
            % Compute difference with respect to 
            err = abs(c_prev - c);
            if err < tol
                break;
            end
            % Evaluate function in middle point
            f_c = fun(c);
            evalcounter = evalcounter + 1;

            % Determine next search interval
            if f_a*f_c < 0
                b = c;
                f_b = f_c;
            else
                a = c;
                f_a = f_c;
            end
            safety_var = safety_var + 1;
            
        end

        fprintf('\nNumber of cycles for bi-section method: %i\n', safety_var);

        % Return output
        x_zero = c;
        feval = f_c;
        info.evalcounter = evalcounter;

    case 2 % Secants

        x_k = b;
        xprev = a;

        f_xprev = f_b;
        f_k = f_a;
        err = 1;

        while abs(err) > tol
            if safety_var == MaxIter
                warning('Solution not found in specified Max Iteration number')
                break;
            end

            x_new = x_k - (x_k - xprev)*f_k./(f_k-f_xprev);

            % Update x and f
            xprev = x_k;
            f_xprev = f_k;

            err = abs(x_new - x_k);

            if err < tol
                break;
            end

            x_k = x_new;       
            f_k = fun(x_k);         
            evalcounter = evalcounter + 1;

            safety_var = safety_var + 1;
        end

        fprintf('\nNumber of cycles for secants method: %i\n', safety_var);

        % Return output
        x_zero = x_k;
        feval = f_k;
        info.evalcounter = evalcounter;

    case 3 % Regula falsi

        x_k = b;
        xprev = a;

        % Initilize variables
        f_xk = f_b;
        f_xprev = f_a;
        err = 1;

        while err > tol

            if safety_var == MaxIter
                break;
            end

            x_new = x_k - (x_k - xprev)*f_xk./(f_xk-f_xprev);

            err = abs(x_k - x_new);

            if err < tol
                break;
            end

            f_xnew = fun(x_new);
            evalcounter = evalcounter + 1;

            if f_xprev*f_xnew < 0
                x_k = x_new;
                f_xk = f_xnew;
            else
                xprev = x_new;
                f_xprev = f_xnew;
            end

            safety_var = safety_var + 1;
        end

        fprintf('\nNumber of cycles for Regula Falsi method: %i \n', safety_var);

        % Return Output
        x_zero = x_new;
        feval = f_xnew;
        info.evalcounter = evalcounter;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCISE 3
function [state, evalcounter, starter_fk] = RKN_integrator(RHS, tspan, state0, order)
%% PROTOTYPE
% [state, evalcounter, starter_fk] = RKN_integrator(RHS, tspan, state0, order)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing RK2 and RK4 schemes to integration vectorial systems
% of ODEs with fixed step-size h. 
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% RHS: [Nx1] function handles specifying the ODEs system. Input variables
%      MUST be specified ast (x, t)
% tspan: [Mx1] vector of times for the integration scheme
% state0: [Kx1] value of the state at initial time t0
% order: [scalar] 2 or 4: selects RK2 or RK4
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% state: [length(tspan), length(state0)] trajectory of the state at tspan times
% evalcounter: [scalar] number of fcn evals
% starter_fk: [length(tspan), length(state0)] function evaluations at
%             points specified by tspan to be used as startup for Multisteps
% -------------------------------------------------------------------------------------------------------------

ifnotstarter = 1;

h_fixed = tspan(2) - tspan(1);
h_vec = h_fixed.*ones(length(tspan), 1);

% Modify last stepsize to solve up to the correct interval
if tspan(end) - tspan(end-1) ~= h_fixed
    h_vec(end) = tspan(end) - tspan(end-1);
end

% If the method is used as starter, fk is evaluated at all point of tspan. Otherwise last cycle is skipper.
if nargout > 2
    ifnotstarter = 0; 
end

% ACHTUNG: RHS MUST BE SPECIFIED WITH VARS AS (x, t)
% NOTE: RHS must be a function with dxdt as output

switch order
    case 2

        % Static memory allocation
        xk = nan(length(tspan)+1, length(state0));
        evalcounter = 0;
        % Initial condition
        xk(1, :) = state0;

        if nargout > 2
            starter_fk = nan(length(tspan), length(state0));
        end

        for idt = 1:length(tspan)

            fk = RHS(xk(idt, :), tspan(idt));
            evalcounter = evalcounter + 1;

            if nargout > 2
                starter_fk(idt, :) = fk;
            end

            % Predictor step
            x_p = xk(idt, :) + h_vec(idt)*fk;
            % Corrector step to get x at k+1 step
            xk(idt + 1, :) = xk(idt, :) + h_vec(idt)* (0.5*fk + 0.5* RHS(x_p, tspan(idt)+h_vec(idt)));
            evalcounter = evalcounter + 1;
        end
        state = xk(1:end-1, :);


    case 4

        %% RK4 Integrator scheme
        %         RK4timer = tic;
        evalcounter = 0;

        time_vec = tspan; %et0:h:etf;
        %         h = tspan(2) - tspan(1); % Fixed step derived from tspan

        alfa = [1/2, 1/2, 1, 1]';

        beta = [1/2, 0, 0, 0;
            0, 1/2, 0, 0;
            0, 0, 1, 0;
            1/6, 1/3, 1/3, 1/6];


        state = nan(length(time_vec), length(state0));

        % Initial condition assignment
        state(1, :) = state0; %x0;

        for timestep = 1:length(time_vec) - ifnotstarter

            x_k = state(timestep, :)';
            t_k = time_vec(timestep);

            % 1st Predictor stage
            fk = RHS(x_k, t_k);
            evalcounter = evalcounter + 1;
            x_P1 = x_k + alfa(1)* h_vec(timestep) * fk;  %f(x_k, t_k)

            if nargout > 2
                if timestep == 1
                    starter_fk = nan(length(tspan), length(fk));
                end

                starter_fk(timestep, :) = fk;
            end

            % 2nd Predictor stage
            fP1 = RHS(x_P1, t_k + alfa(1)*h_vec(timestep));
            evalcounter = evalcounter + 1;
            x_P2 = x_k + alfa(2)*h_vec(timestep) *fP1;  %f(x_P1, t_k + alfa(1)*h)

            % 3rd Predictor stage
            fP2 = RHS(x_P2,t_k + alfa(2)*h_vec(timestep));
            evalcounter = evalcounter + 1;
            x_P3 = x_k + alfa(3)*h_vec(timestep) *fP2; %f(x_P2, t_k + alfa(2)*h)

            % 4th Corrector stage
            fP3 = RHS(x_P3, t_k + alfa(3)*h_vec(timestep));
            evalcounter = evalcounter + 1;

            x_next = x_k + alfa(4).*h_vec(timestep) .* (beta(4, 1)*fk  + beta(4, 2)* fP1 + beta(4, 3)* fP2 + beta(4, 4)* fP3);

            state(timestep + 1, :) = x_next;

        end

        if ifnotstarter == 0
            state = state(1:end-1, :);
        end

        %         toc(RK4timer)

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCISE 5
function [xk_tf] = RK_evalfinalstate(F, x0, h, t0, tf)
%% PROTOTYPE
% [xk_tf] = RK_evalfinalstate(F, x0, h, t0, tf)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes state x_k at final time tf by using the forward operator
% specified by F, starting from state x0
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% F: [function handle] specified with input variable h only
% x0: [Mx1] value of state at initial time
% h: [scalar] stepsize for the integration
% t0, tf: [scalars] initial and final times
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------------------------------------------

% Evaluate how many nsteps are required to cover tspan
h = abs(h);
nstep = floor( (tf-t0)/h ); 

if h*nstep > tf
    warning('tf has been exceeded')
end

% Compute final state of integration
if h*nstep == tf
    xk_tf = (F(h)^nstep)*x0;
else
    
    xk_tf = F(tf - h*(nstep))*(F(h)^nstep)*x0;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCISE 7
function [x, evalcounter] = IEX4_integrator(RHS, tspan, x0)
%% PROTOTYPE
% [x, evalcounter] = IEX4_integrator(RHS, tspan, x0)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Integrator implementing IXE4 scheme for vectorial ODEs at fixed step size
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% RHS: [Nx1] function handle specifiyng the dynamics. MUST be defined with
%       variables (t, x) in this order
% tspan: [Mx1] vector of times for the integration
% x0: [Kx1] vector of states at initial time tspan(1)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% x: [length(tspan), K] trajectory of the state at each time point
% evalcounter: [scalar] counter of function evalutions
% -------------------------------------------------------------------------------------------------------------

%% Function code

IEX4_timer = tic;

% Integrator parameters and counters
alfa = [-1/6, 4, -27/2, 32/3]';
h_fixed = tspan(2) - tspan(1);
h = [1, 1/2, 1/3, 1/4] * h_fixed;

options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-12, 'UseParallel', false, 'Algorithm', 'levenberg-marquardt');

% Static memory allocation
x = nan(length(tspan), length(x0));
x(1, :) = x0;

evalcounter = 0;

for timestep = 1:length(tspan)-1

    x_k = x(timestep, :)';
    t_k = tspan(timestep);

    % 1st Predictor stage
    % Note: stages evaluations are independent: can be parallelized

    % Implicit non-linear equation to solve
%     f = @(xP1) xP1 - x_k - h(1) * RHS(t_k + h(1), xP1);
    [xP1, ~, ~, output] = fsolve(@(xP1) xP1 - x_k - h(1) * RHS(t_k + h(1), xP1), x_k, options);  
    evalcounter = evalcounter + output.funcCount;

    % 2nd Predictor stage
%     f = @(xP2_mid) xP2_mid - x_k - h(2) * RHS(t_k + h(2), xP2_mid);
    [xP2_a, ~, ~, output] = fsolve(@(xP2_a) xP2_a - x_k - h(2) * RHS(t_k + h(2), xP2_a), x_k, options);
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP2) xP2 - xP2_a - h(2) * RHS(t_k + h(1), xP2); % Nbp_sys(t_k + h(1), xP2, bodies, Rframe, Oframe);
    [xP2, ~, ~, output] = fsolve(@(xP2) xP2 - xP2_a - h(2) * RHS(t_k + h(1), xP2), xP2_a, options);  
    evalcounter = evalcounter + output.funcCount;

    % 3th Predictor stage

%     f = @(xP3_a) xP3_a - x_k - h(3) * RHS(t_k + h(3), xP3_a);
    [xP3_a, ~, ~, output] = fsolve(@(xP3_a) xP3_a - x_k - h(3) * RHS(t_k + h(3), xP3_a), x_k, options); 
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP3_b) xP3_b - xP3_a - h(3) * RHS(t_k + h(3), xP3_b);
    [xP3_b, ~, ~, output] = fsolve(@(xP3_b) xP3_b - xP3_a - h(3) * RHS(t_k + h(3), xP3_b), xP3_a, options);
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP3) xP3 - xP3_b - h(3) * RHS(t_k + h(1), xP3); % Nbp_sys(t_k + h(1), xP3, bodies, Rframe, Oframe);
    [xP3, ~, ~, output] = fsolve(@(xP3) xP3 - xP3_b - h(3) * RHS(t_k + h(1), xP3), xP3_b, options);
    evalcounter = evalcounter + output.funcCount;

    % 4th Predictor stage

%     f = @(xP4_a) xP4_a - x_k - h(4) * RHS(t_k + h(4), xP4_a);
    [xP4_a, ~, ~, output] = fsolve(@(xP4_a) xP4_a - x_k - h(4) * RHS(t_k + h(4), xP4_a), x_k, options); 
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP4_b) xP4_b - xP3_a - h(4) * RHS(t_k + h(4), xP4_b);
    [xP4_b, ~, ~, output] = fsolve(@(xP4_b) xP4_b - xP4_a - h(4) * RHS(t_k + h(4), xP4_b), xP4_a, options);
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP4_c) xP4_c - xP4_b - h(4) * RHS(t_k + h(4), xP4_c);
    [xP4_c, ~, ~, output] = fsolve(@(xP4_c) xP4_c - xP4_b - h(4) * RHS(t_k + h(4), xP4_c), xP4_b, options);
    evalcounter = evalcounter + output.funcCount;

%     f = @(xP4) xP4 - xP4_c - h(4) * RHS(t_k + h(1), xP4); % Nbp_sys(t_k + h(1), xP4, bodies, Rframe, Oframe);
    [xP4, ~, ~, output] = fsolve(@(xP4) xP4 - xP4_c - h(4) * RHS(t_k + h(1), xP4), xP4_c, options);
    evalcounter = evalcounter + output.funcCount;

    % 5th Corrector stage
    x_next = alfa(1)*xP1 + alfa(2)*xP2 + alfa(3)*xP3 + alfa(4)*xP4;

    x(timestep + 1, :) = x_next;
end

toc(IEX4_timer)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 8
function [time_vec, x, evalcounter] = Multistep_integrator(RHS, tspan, x0, class, ~)
%% PROTOTYPE
% [time_vec, x, evalcounter] = Multistep_integrator(RHS, tspan, x0, class, ~)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function implementing multistep methods: AB3, ABM3, AM3, BDF3. Not yet
% implemented: AM4, AB4, ABM4. RK4 is used as startup method.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% RHS: [Nx1] fcn handle specifying the dynamics. MUST be defined with
%       variables (x, t) in this order
% tspan: [Mx1] vector of times for the integration
% x0: [Kx1] vector of states at initial time tspan(1)
% class: [scalar] specifying the scheme to be used. 1) AB, 2) AM, 3) ABM, 4) BDF
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% time_vec: [Mx1] equal to tspan: vector of time points
% x: [MxK] trajectory of the state at each time point
% evalcounter: [scalar] number of fcn evals
% -------------------------------------------------------------------------------------------------------------

fsolve_tol = 1e-12;
order = 3;

%% Function code
% ACHTUNG: RHS MUST BE SPECIFIED WITH VARS AS (x, t)

% Scheme selection
if strcmpi('AB', class)
    class = 0;
elseif strcmpi('AM', class)
    class = 1;
elseif strcmpi('ABM', class)
    class = 2;
elseif strcmpi('BDF', class)
    class = 3;
else
    err('Specified method name not found or available')
end


% Schemes parameters and counters
h = tspan(2) - tspan(1);
starting_idt = 1;
starter_counter = 0;
evalcounter = 0;

% Static memory allocation
xk = nan(length(tspan), length(x0)); % Allocation memory for solution

fk = nan(length(tspan), length(RHS(x0, tspan(1))));
evalcounter = evalcounter + 1;

xk(1, :) = x0;


if class == 0 || class == 1 || class == 2 % AM/AB/ABM

    switch order
        case 3  
            [x_starter, starter_counter, starter_fk] = RKN_integrator(RHS, tspan(1:3), x0, 4);
            xk(1:3, :) = x_starter;
            fk(1:3, :) = starter_fk;

            if class == 1
                starting_idt = 2;
            else
                starting_idt = 3;
            end
            
        case 4  
            [x_starter, starter_counter, starter_fk] = RKN_integrator(RHS, tspan(1:4), x0, 4);
            xk(1:3, :) = x_starter;
            fk(1:3, :) = starter_fk;
            starting_idt = 4;
    end

elseif class == 3 % BDF

    switch order
        case 3 % BDF3
            [x_starter, starter_counter, starter_fk] = RKN_integrator(RHS, tspan(1:3), x0, 4);
            xk(1:3, :) = x_starter;
            fk(1:3, :) = starter_fk;
            
            starting_idt = 3;
    end
end

evalcounter = evalcounter + starter_counter;

switch class

    case 0 % 'AB'

        for idt = starting_idt:length(tspan)-1

            tk = tspan(idt);          
            fk(idt, :) = RHS(xk(idt, :), tk);
            evalcounter = evalcounter + 1;

            if order == 3
                % AB3
                xk(idt+1, :) = xk(idt, :) + 1/12 * h*(23*fk(idt, :) - 16*fk(idt-1, :) + 5*fk(idt-2, :));

            elseif order == 4
                % AB4
                xk(idt+1, :) = xk(idt, :) + 1/24 * h*(55*fk(idt, :) - 59*fk(idt-1, :) + 37*fk(idt-2) - 9*fk(idt-3, :));

            end
        end

    case 1 %'AM'

        opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', fsolve_tol);

        for idt = starting_idt:length(tspan)-1

            tk = tspan(idt);

            fk(idt, :) = RHS(xk(idt, :), tk);
            evalcounter = evalcounter + 1;

            if order == 3
                % AM3
                [xk(idt+1, :), ~, ~, output] = fsolve(@(x_next) x_next - xk(idt, :) - (1/12 * h* ( 5*(RHS(x_next, tk+h))' + 8*fk(idt, :) - fk(idt-1, :) ) ) , xk(idt, :), opts);
                evalcounter = evalcounter + output.funcCount;

            elseif order == 4
                % AM4
                [xk(idt+1, :), ~, ~, output] = fsolve(@(x_next) x_next - xk(idt, :) - (1/24 * h*(9*(RHS(x_next, tk+h))' + 19*fk(idt, :) - 5*fk(idt-1, :) + fk(idt-2, :))) , xk(idt, :), opts);
                evalcounter = evalcounter + output.funcCount;
            end
        end

    case 2 %'ABM'
        for idt = starting_idt:length(tspan)-1

            tk = tspan(idt);
            % COMPUTE idt-1 and idt-2 with RK3 or 4
            fk(idt, :) = RHS(xk(idt, :), tk);
            evalcounter = evalcounter + 1;

            if order == 3
                % ABM3
                % Predictor stage - AB3
                xp = xk(idt, :) + 1/12 * h*(23*fk(idt, :) - 16*fk(idt-1, :) + 5*fk(idt-2, :));
                % Corrector stage - AM3
                xk(idt+1, :) = xk(idt, :) + 1/12 * h*(5*(RHS(xp, tk+h))' + 8*fk(idt, :) - fk(idt-1, :));
                evalcounter = evalcounter + 1;
            elseif order == 4
                % ABM4
                % Predictor stage - AB4
                xp = xk(idt, :) + 1/24 * h*(55*fk(idt, :) - 59*fk(idt-1, :) + 37*fk(idt-2) - 9*fk(idt-3, :));
                % Corrector stage - AM4
                xk(idt+1, :) = xk(idt, :) + 1/24 * h*(9*(RHS(xp, tk+h))' + 19*fk(idt, :) - 5*fk(idt-1, :) + fk(idt-2, :));
                evalcounter = evalcounter + 1;
            end
        end

    case 3 %'BDF'

        % BDF3: xk(idt+1, :) = 18/11 * xk(idt, :) - 9/11 * xk(idt-1, :) + 2/11 * xk(idt-2, :) + 6/11* h * RHS(x_next, tk+h);

        opts = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', fsolve_tol);

        for idt = starting_idt:length(tspan)-1

            tk = tspan(idt);

            [xk(idt+1, :), ~, ~, output] = fsolve(@(x_next) x_next - (18/11 * xk(idt, :) - 9/11 * xk(idt-1, :) + 2/11 * xk(idt-2, :) + 6/11* h * (RHS(x_next, tk+h))') , xk(idt, :), opts);
            evalcounter = evalcounter + output.funcCount;
        end

end

time_vec = tspan;
x = xk;
end
