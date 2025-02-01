function [R_313] = EA313toDCM(phi, theta, psi, unit)

%% PROTOTYPE
% [R_313] = EA_313toDCM(phi, theta, psi, unit)
% ---------------------------------------------------------------------
%% DESCRIPTION
% Remark: the multiplication order is inverted wrt physical order of rotation
% R_313(phi, theta, psi) = R_3(psi)*R_1(theta)*R_3(phi)
% ---------------------------------------------------------------------
%% INPUT
% phi: [rad] 1st rotation angle about Z axis
% theta: [rad] 2nd rotation angle about X axis
% psi: [rad] 3rd rotation angle about Z axis
% unit: [string] Specifies the unit of measure. DEFAULT: rad
% ---------------------------------------------------------------------
%% OUTPUT
% ---------------------------------------------------------------------
%% CHANGELOG
% V1 coded: basic conversion, no check for singularities, 18/10/2021
% V1 discarded: not working properly, 21/10/2021
% V2 coded: basic conversion w/ choice of unit of measure, 21/10/2021
% ---------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano

%% Next upgrades
% 1) Check for singularities
 
if nargin < 4

unit = "rad";

end

R_313 = Rot3(psi, unit) * Rot1(theta, unit) * Rot3(phi, unit);

end



