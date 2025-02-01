function [domain] = nonlinspace(LowerB, UpperB, exponent, N_points, twosided)
%% PROTOTYPE
% [domain] = nonlinspace(LowerB, UpperB, exponent, N_points, twosided)
% -------------------------------------------------------------------------
%% INPUT
% LowerB: lower bound of the domain
% UpperB: upper bound of the domain
% exponent: control parameter to change the density of the distribution
% N_points: number of points for the domain
% twosided: 0 or 1 --> 0: one sided, 1: two sided
% -------------------------------------------------------------------------
%% OUTPUT
% domain: non linearly spaced array
% -------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano
% -------------------------------------------------------------------------
%% CHANGELOG
% V1: generates one sided and two sided non linearly spaced array, 14/12/2021
% -------------------------------------------------------------------------

%% Function code

alfa = exponent; % decides where and how mush the point density is increa
csi = linspace(LowerB^(1/alfa), UpperB^(1/alfa), N_points); % uniformly spaced mesh to be passed to the mappng function
x = csi.^alfa; % mapping function

if twosided == 1
    x = [-flip(x) x]; % Se voglio che il mesh sia simmetrico e denso ai lati
    y = x - UpperB; % translate the domain to [0, UpperB]
    domain = y;
else
    domain = x;
end

end
