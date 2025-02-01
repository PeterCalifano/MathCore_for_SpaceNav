function o_dDCM = Quat2DCM(i_dQuatRot, i_bIS_VSRPplus) %#codegen
arguments
    i_dQuatRot     (4,1) double
    i_bIS_VSRPplus (1,1) logical = true
end
%% PROTOTYPE
% o_dDCM = Quat2DCM(i_dQuatRot, i_bIS_VSRPplus) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting attitude quaternion to DCM. Conversion occurs according to VSRP+ convention
% (right-handed) Set i_bIS_VSRPplus = false if SVRP+ convention with scalar first is being used. 
% (SV) Scalar first, Vector last
% (P) Passive 
% (R) Successive coordinate transformations have the unmodified quaternion chain on the Right side of
%     the triple product.
% (plus) Right-Handed Rule for the imaginary numbers i, j, k. (aka Hamilton)
% REFERENCE:
% 1) Section 2.9.3, page 45, Fundamentals of Spacecraft Attitude 
% Determination and Control, Markley, Crassidis, 2014
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dQuatRot:       [4, 1] Input quaternion to convert to DCM
% i_bIS_VSRPplus:   [1]  Boolean flag indicating the convention of the
%                       quaternion. 1: VSRPplus, 0: SVRPplus (as MATLAB)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dDCM: [3, 3] Output rotation matrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 01-09-2023    Pietro Califano     Coded from reference. Validated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% DCM initialization
o_dDCM = coder.nullcopy(zeros(3,3));

% Get quaternion components depending on quaternion order
if i_bIS_VSRPplus == false
    % SVRPplus convention (true)
    qs  = i_dQuatRot(1);
    qv1 = i_dQuatRot(2);
    qv2 = i_dQuatRot(3);
    qv3 = i_dQuatRot(4);
else
    % VSRPplus convention (false)
    qv1 = i_dQuatRot(1);
    qv2 = i_dQuatRot(2);
    qv3 = i_dQuatRot(3);
    qs  = i_dQuatRot(4);
end

% Convert to DCM 
o_dDCM(1,1) = qs^2 + qv1^2 - qv2^2 - qv3^2;
o_dDCM(2,1) = 2*(qv1*qv2 - qv3*qs);
o_dDCM(3,1) = 2*(qv1*qv3 + qv2*qs);

o_dDCM(1,2) = 2*(qv1*qv2 + qv3*qs);
o_dDCM(2,2) = qs^2 - qv1^2 + qv2^2 - qv3^2;
o_dDCM(3,2) = 2*(qv2*qv3 - qv1*qs);

o_dDCM(1,3) = 2*(qv1*qv3 - qv2*qs);
o_dDCM(2,3) = 2*(qv2*qv3 + qv1*qs);
o_dDCM(3,3) = qs^2 - qv1^2 - qv2^2 + qv3^2;

end
