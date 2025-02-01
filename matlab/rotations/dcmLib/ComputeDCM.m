function o_dRotDCM = ComputeDCM(i_ui8AxisIDseq, i_dRotAngleSeq) %#codegen
%% PROTOTYPE
% o_dRotDCM = ComputeDCM(i_ui8AxisIDseq, i_dRotAngleSeq)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function computing a rotation matrix by concatenation of the selected
% single axis rotation sequence of N entries (arbitrary number).
% Angles must be in radians. Rotation order is "as written": first in the
% sequence is the last rotation to be perfomed.
% NOTE: Rotation matrix assumed Right-Hand Rule, positive clockwise in the
%       direction of the axis.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8AxisIDseq: [1, N] Array sequence of axis IDs (1, 2 or 3)
% i_dRotAngleSeq: [1, N] Array sequence of corresponding rotation angles
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dRotDCM: [3, 3] Output rotation matrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 27-11-2023    Pietro Califano     Function assembled from validated code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Get number of rotations
Nr = length(i_dRotAngleSeq);

% Input checks
assert(Nr == length(i_ui8AxisIDseq), 'Input arrays length do not match.')
assert( sum(i_ui8AxisIDseq >= 1 & i_ui8AxisIDseq <= 3) == Nr , 'Invalid axis selected. Input must be: 1, 2 or 3.');

% Initialize DCM with first rotation (left-most)
o_dRotDCM = SingleAxisDCM(i_ui8AxisIDseq(end), i_dRotAngleSeq(end));

for idR = (Nr-1):-1:1
    % Concatenate with with next rotation matrix
    o_dRotDCM = SingleAxisDCM(i_ui8AxisIDseq(idR), i_dRotAngleSeq(idR)) * o_dRotDCM;
end


%% LOCAL FUNCTIONS
    function o_dSingleRotDCM = SingleAxisDCM(i_ui8AxisID, i_dRotAngle) %#codegen
        %% PROTOTYPE
        % o_dSingleRotDCM = SingleAxisDCM(i_ui8AxisID, i_dRotAngle)
        % -------------------------------------------------------------------------------------------------------------
        %% DESCRIPTION
        % Function computing the rotation matrix about the selected axis (1,2,3) =
        % (X,Y,Z) of the angle i_dRotAngle in RADIANS. Zero is returned if invalid
        % axis input is selected and a warning is issued.
        % NOTE: Rotation matrix assumed Right-Hand Rule, positive clockwise in the
        %       direction of the axis.
        % -------------------------------------------------------------------------------------------------------------
        %% INPUT
        % i_ui8AxisID: [1] Axis ID (1,2,3) = (X,Y,Z)
        % i_dRotAngle: [1] Rotation angle in RADIANS
        % -------------------------------------------------------------------------------------------------------------
        %% OUTPUT
        % o_dSingleRotDCM: [3, 3] DCM rotating around "AxisID" of "RotAngle" [rad]
        % -------------------------------------------------------------------------------------------------------------
        %% CHANGELOG
        % 27-11-2023    Pietro Califano     Function assembled from
        %                                   previous validated codes.
        % -------------------------------------------------------------------------------------------------------------
        %% DEPENDENCIES
        % [-]
        % -------------------------------------------------------------------------------------------------------------
        %% Function code
        switch i_ui8AxisID
            case 1
                % Axis 1: X
                o_dSingleRotDCM = [1, 0, 0; ...
                    0, cos(i_dRotAngle), sin(i_dRotAngle); ...
                    0, -sin(i_dRotAngle), cos(i_dRotAngle)];
            case 2
                % Axis 2: Y
                o_dSingleRotDCM = [cos(i_dRotAngle), 0, -sin(i_dRotAngle);...
                    0, 1, 0;...
                    sin(i_dRotAngle), 0, cos(i_dRotAngle)];
            case 3
                % Axis 3: Z
                o_dSingleRotDCM = [cos(i_dRotAngle), sin(i_dRotAngle), 0;...
                    -sin(i_dRotAngle), cos(i_dRotAngle), 0;...
                    0, 0, 1];
            otherwise
                disp('Invalid rotation axis ID. Output set to zero!')
                o_dSingleRotDCM = zeros(3, 3);
        end
    end

end
