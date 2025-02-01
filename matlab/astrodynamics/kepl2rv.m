function xCart = kepl2rv(xKepl, mu, unit) %#codegen
%% PROTOTYPE
% xCart = kepl2rv(xKepl, mu, unit)
% --------------------------------------------------------------------------
%% INPUT
% xKepl = [6x1] Vector of Keplerian elements 
%   SMA: [1] Semi Major axis of the orbit
%   ECC: [1] Eccentricity magnitude
%   incl: [1] Inclination of orbital plane wrt XY plane (dihedral angle)
%   RAAN: [1] Angle between the Node Line and the X axis direction 
%   omega: [1] Angle pinpointing eccentricity vector from Node Line
%   TA: [1] True anomaly theta
% unit: [string] "deg" if input angles in deg. (default: rad)
% --------------------------------------------------------------------------
%% OUTPUT
% xCart: [6x1] Cartesian Coordinates from Osculating elements xKepl
% --------------------------------------------------------------------------
%% CHANGELOG
% V1 coded, not tested yet: basic conversion, no check for singularities, 18/10/2021. 
% V1 not working correctly, 22/10/2021
% V1 debugged, now working: error in the building blocks of DCM, 23/10/2021
% 12-07-2023    Pietro Califano     Re-worked and validated
% --------------------------------------------------------------------------


%% Input pre-processing
SMA = xKepl(1);
ECC = xKepl(2);
incl = xKepl(3);
RAAN = xKepl(4);
omega = xKepl(5);
TA = xKepl(6);

if nargin < 3
    unit = 'rad';
end

if unit == "deg"
    TA = deg2rad(TA);
    incl = deg2rad(incl);
    RAAN = deg2rad(RAAN);
    omega = deg2rad(omega);
end

%% Conversion [SMA, ECCTA] --> (R, V) in perifocal frame
% The semi-major axis a, the eccentricity e and the true anomaly defines
% the shape (so the energy) of the orbit; true anomaly identifies the
% position at the specific time t of the satellite along the orbit. 
% From the first two, h can be computed. Then, r and v in perifocal coords.
% NOTE: the third components is 0 by definition.

% Norm of the specific orbital ang. mom.
h_norm = sqrt(mu*SMA*(1 - ECC^2)); 

r_norm = (h_norm.^2/mu) .* (1/(1 + ECC.*cos(TA)));

r_perifocal = r_norm .* [cos(TA); sin(TA); 0]; % 3x1
v_perifocal = [-(mu./h_norm).*sin(TA); (mu./h_norm ).*(ECC + cos(TA)); 0]; % 3x1

%% DCM Perifocal (e, p, h) --> Cartesian Inertial (I, J, K)
% Having r and v in perifocal and the Euler angles (313 seq) that defines the orbital 
% plane wrt Cartesian Inertial, the transformation is done by means of the
% corresponding DCM (ZXZ). Inclination rotation matches the XY plane.

% Order of rotation: 1st RAAN, 2nd inclination, 3rd omega
ECEI2perifocal = Rot3(omega, "rad") * Rot1(incl, "rad") * Rot3(RAAN, "rad");

xCart = zeros(6, 1);
xCart(1:3) = ECEI2perifocal' * r_perifocal; 
xCart(4:6) = ECEI2perifocal' * v_perifocal;

%% LOCAL functions
function [Rot_Mat1] = Rot1(rot_angle, unit)

%% PROTOTYPE
% [Rot_Mat1] = Rot1(rot_angle)
% ----------------------------------------------------------------------------------------
%% DESCRIPTION
% The following code build the matrix representing a rotation of given angle around axis 1 
% ----------------------------------------------------------------------------------------
%% INPUT
% rot_angle: [scalar] angle of rotation in rad or deg
% unit: [string] string to specify the unit of measure of the angle: can be
%        rad or deg. DEFAULT: rad.
% ----------------------------------------------------------------------------------------
%% OUTPUT
% Rot_mat1 : [3x3] rotation matrix
% ----------------------------------------------------------------------------------------
%% CHANGELOG
% V1 coded: complete, 21/10/2021
% ----------------------------------------------------------------------------------------
%% CONTRIBUTORS
% Pietro Califano

if nargin == 1
    unit = "rad";
end

if unit == "rad"

    Rot_Mat1 = [1, 0, 0; ...
        0, cos(rot_angle), sin(rot_angle); ...
        0, -sin(rot_angle), cos(rot_angle)];

elseif unit == "deg"

    Rot_Mat1 = [1, 0, 0; ...
        0, cosd(rot_angle), sind(rot_angle); ...
        0, -sind(rot_angle), cosd(rot_angle)];
end

end

    function [Rot_Mat3] = Rot3(rot_angle, unit)

        %% PROTOTYPE
        % [Rot_Mat3] = Rot3(rot_angle)
        % ----------------------------------------------------------------------------------------
        %% DESCRIPTION
        % The following code build the matrix representing a rotation of given angle around axis 3
        % ----------------------------------------------------------------------------------------
        %% INPUT
        % rot_angle: [scalar] angle of rotation in rad or deg
        % unit: [string] string to specify the unit of measure of the angle: can be
        %        rad or deg. DEFAULT: rad.
        % ----------------------------------------------------------------------------------------
        %% OUTPUT
        % Rot_mat3: [3x3] rotation matrix
        % ----------------------------------------------------------------------------------------
        %% CHANGELOG
        % V1 coded: complete, 21/10/2021
        % ----------------------------------------------------------------------------------------
        %% CONTRIBUTORS
        % Pietro Califano

        if nargin == 1
            unit = "rad";
        end

        if unit == "rad"

            Rot_Mat3 = [cos(rot_angle), sin(rot_angle), 0;...
                -sin(rot_angle), cos(rot_angle), 0;...
                0, 0, 1];

        elseif unit == "deg"

            Rot_Mat3 = [cosd(rot_angle), sind(rot_angle), 0;...
                -sind(rot_angle), cosd(rot_angle), 0;...
                0, 0, 1];

        end
    end

end