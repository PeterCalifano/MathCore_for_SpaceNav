function dDCM = Quat2DCM(dQuatRot, bIS_VSRPplus) %#codegen
arguments
    dQuatRot     (4,1) double
    bIS_VSRPplus (1,1) logical = true
end
%% PROTOTYPE
% dDCM = Quat2DCM(dQuatRot, bIS_VSRPplus) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting attitude quaternion to DCM. Conversion occurs according to VSRP+ convention
% (right-handed) Set bIS_VSRPplus = false if SVRP+ convention with scalar first is being used. 
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
% dQuatRot:       [4, 1] Input quaternion to convert to DCM
% bIS_VSRPplus:   [1]  Boolean flag indicating the convention of the
%                       quaternion. 1: VSRPplus, 0: SVRPplus (as MATLAB)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dDCM: [3, 3] Output rotation matrix
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
dDCM = coder.nullcopy(zeros(3,3));

% Get quaternion components depending on quaternion order
if bIS_VSRPplus == false
    % SVRPplus convention (true)
    qs  = dQuatRot(1);
    qv1 = dQuatRot(2);
    qv2 = dQuatRot(3);
    qv3 = dQuatRot(4);
else
    % VSRPplus convention (false)
    qv1 = dQuatRot(1);
    qv2 = dQuatRot(2);
    qv3 = dQuatRot(3);
    qs  = dQuatRot(4);
end

% Convert to DCM 
dDCM(1,1) = qs^2 + qv1^2 - qv2^2 - qv3^2;
dDCM(2,1) = 2*(qv1*qv2 - qv3*qs);
dDCM(3,1) = 2*(qv1*qv3 + qv2*qs);

dDCM(1,2) = 2*(qv1*qv2 + qv3*qs);
dDCM(2,2) = qs^2 - qv1^2 + qv2^2 - qv3^2;
dDCM(3,2) = 2*(qv2*qv3 - qv1*qs);

dDCM(1,3) = 2*(qv1*qv3 - qv2*qs);
dDCM(2,3) = 2*(qv2*qv3 + qv1*qs);
dDCM(3,3) = qs^2 - qv1^2 - qv2^2 + qv3^2;

end
