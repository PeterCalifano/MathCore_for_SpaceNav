function o_dQuatRot = DCM2quat(i_dDCM, i_bIS_VSRPplus) %#codegen
%% PROTOTYPE
% o_dQuatRot = DCM2quat(i_dDCM, i_bIS_VSRPplus)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting a DCM to Attitude quaternion. Conversion occurs according to VSRP+ convention
% (right-handed) Set i_bIS_VSRPplus = false if SVRP+ convention with scalar first is being used. 
% The conversion is made numerically optimal. Only supports one conversion per call.
% (SV) Scalar first, Vector last
% (P) Passive 
% (R) Successive coordinate transformations have the unmodified quaternion chain on the Right side of
%     the triple product.
% (plus) Right-Handed Rule for the imaginary numbers i, j, k. (aka Hamilton)
% REFERENCE:
% 1) Section 2.9.3, page 48, Fundamentals of Spacecraft Attitude
% Determination and Control, Markley, Crassidis, 2014
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% o_dDCM:         [3, 3]  Input rotation matrix to convert
% i_bIS_VSRPplus: [1]  Boolean flag indicating the convention of the
%                      quaternion. 1: VSRPplus, 0: SVRPplus (as MATLAB)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dQuatRot:     [4, 1]  Output quaternion
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-12-2023    Pietro Califano     Coded from reference. Validated against
%                                   MATLAB functions.
% 27-04-2024    Pietro Califano     Modified convention of quaternion
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% INPUT ASSERT CHECKS
assert( all(size(i_dDCM) == [3,3], 'all'), 'ERROR: DCM must be [3x3x1]')

% Quaternion output initialization
o_dQuatRot = coder.nullcopy(zeros(4, 1));

% Compute trace of DCM
trDCM = trace(i_dDCM);

[~, idMax] = max([i_dDCM(1,1), i_dDCM(2,2), i_dDCM(3,3), trDCM]);

% Adjust computation case for numerical accuracy to convert DCM to NON
% UNITARY QUATERNION
switch idMax
    case 4
        % Scalar part is max
        o_dQuatRot = [i_dDCM(2, 3) - i_dDCM(3, 2);
            i_dDCM(3, 1) - i_dDCM(1, 3);
            i_dDCM(1, 2) - i_dDCM(2, 1);
            1 + trDCM];
    case 3
        % Third vector component is max
        o_dQuatRot = [i_dDCM(3, 1) + i_dDCM(1, 3);
            i_dDCM(3, 2) + i_dDCM(2, 3);
            1 + 2*i_dDCM(3, 3) - trDCM;
            i_dDCM(1, 2) - i_dDCM(2, 1)];

    case 2
        % Second vector component is max
        o_dQuatRot = [i_dDCM(2, 1) + i_dDCM(1, 2);
            1 + 2*i_dDCM(2, 2) - trDCM;
            i_dDCM(2, 3) + i_dDCM(3, 2);
            i_dDCM(3, 1) - i_dDCM(1, 3)];

    case 1
        % First vector component is max
        o_dQuatRot = [1 + 2*i_dDCM(1, 1) - trDCM;
            i_dDCM(1, 2) + i_dDCM(2, 1);
            i_dDCM(1, 3) + i_dDCM(3, 1);
            i_dDCM(2, 3) - i_dDCM(3, 2)];
end

% Normalize quaternion
o_dQuatRot = o_dQuatRot./norm(o_dQuatRot);

% Assign quaternion components depending on quaternion convention
if i_bIS_VSRPplus == false
    % SVRPplus convention scalar first (SWAP COMPONENTS)
    o_dQuatRot = [o_dQuatRot(4); o_dQuatRot(1:3)];

elseif i_bIS_VSRPplus == true
    % VSRPplus convention (false). DO NOTHING: DEFAULT
    %     o_dQuatRot(1) = qv1;
    %     o_dQuatRot(2) = qv2;
    %     o_dQuatRot(3) = qv3;
    %     o_dQuatRot(4) = qs ;
end

end


