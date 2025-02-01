function xKepl = rv2kepl(xCart, mu, unit) %#codegen
%% PROTOTYPE
% xKepl = rv2oe(xCart, mu, unit)
%   with xKepl = [SMA; ECC; incl; RAAN; omega; TA];
% --------------------------------------------------------------------------
%% DESCRIPTION
% Executes the conversion from Position-Velocity State to Keplerian
% elements. Singularities are managed setting RAAN = 0 if i = 0 (equatorial
% orbit) and argument_of_pericentre = 0 if e = 0 (circular orbit);
% All angles are given as output in [-180, 180] degree interval.
% --------------------------------------------------------------------------
%% INPUT
% xCart: [6x1] State vector in Cartesian coordinates
% mu: Gravitational Parameter of the Main attractor
% --------------------------------------------------------------------------
%% OUTPUT
% xKepl = [6x1] Vector of Keplerian elements 
%   SMA: [1] Semi Major axis of the orbit [unit]
%   ECC: [1] Eccentricity magnitude [unit]
%   incl: [1] Inclination of orbital plane wrt XY plane (dihedral angle) [unit]
%   RAAN: [1] Angle between the Node Line and the X axis direction [unit]
%   omega: [1] Angle pinpointing eccentricity vector from Node Line [unit]
%   TA: [1] True anomaly theta [unit]
% ------------------------------------------------------------- -------------
%% CHANGELOG
% V1 coded, not tested yet: basic conversion, no check for singularities, 18/10/2021.
% V1 tested: 19/10/2021
% V1 not working correctly, 22/10/2021
% V1 debug correction: wrong inverse trigonometric functions, 22/10/2021
% V1 debugged: properly working now, 23/10/2021
% V2: added code to manage singularity cases, 23/11/2021
% V3: debugged omega calculation, 08/12/2021
% 12-07-2023    Pietro Califano     New version, validated and checks for
%                                   Hyperbolic orbits
% --------------------------------------------------------------------------


%% Function code
% Initialize variables
SMA = -1;
ECC = -1;
incl = -1;
RAAN = -1;
omega = -1;
TA = -1;


r = xCart(1:3);
v = xCart(4:6);
r_norm = norm(r); % [LU]
v_norm = norm(v); % [LU/TU]

% Specific orbital angular momentum
h_vector = cross(r, v);
h_norm = norm(h_vector);

% Vr Radial velocity
Vr = dot(r, v)/r_norm;

% 1) Inclination
incl = acos(h_vector(3)./h_norm);

% Node line unit vetor
N = cross([0, 0, 1], h_vector);

% 2) Right Ascension of Ascending Node
if incl <= 1e-12
    % Manage inclination singularity 
    % Set node line coincident with the X-axis
    RAAN = 0; 
    N_line = [1, 0, 0]; 
else
    N_line = N./norm(N);

    if N_line(2) >= 0
        RAAN = acos(N_line(1));
    elseif  N_line(2) < 0
        RAAN = - acos(N_line(1));
    end
end

% Eccentricity vector 
e_vec = (1/mu)*((v_norm.^2 - mu./r_norm) * r - r_norm .* Vr .* v);
% e_vecCheck = cross(v, h_vector)/mu - r/r_norm;

% 3) Eccentricity
ECC = norm(e_vec);
e_unit = e_vec./ECC;


% 4) Argument of Pericentre
if ECC <= 1e-12
    omega = 0;
else
    if  e_vec(3) > 0
        omega = acos(dot(N_line, e_unit));
    elseif  e_vec(3) <= 0
        omega = - acos(dot(N_line, e_unit));
    end
end

% 5) Semi major axis
SMA = (h_norm.^2)/(mu*(1 - ECC.^2));

% 6) True anomaly 
e_dot_r = dot(e_unit, r./r_norm);

if Vr >= 0
    TA = acos(e_dot_r);
    % Assert for hyperbolic orbits

    if ECC > 1
        TAinf = pi - acos(1/ECC);
        assert((abs(TA) - abs(TAinf) < 1e-5), 'TA resulted greater than TA limit!')
    end

else
    TA = - acos(e_dot_r);
    % Assert for hyperbolic orbits

    if ECC > 1
        TAinf = pi - acos(1/ECC);
        assert((abs(TA) - abs(TAinf) < 1e-5), 'TA resulted greater than TA limit!')
    end

end


% Stack and convert from rad to deg if needed
xKepl = [SMA; ECC; incl; RAAN; omega; TA];

if nargin < 4
    unit = 'rad';
end

if strcmpi(unit, 'deg')
    xKepl(3:6) = rad2deg(xKepl(3:6));
end

end










