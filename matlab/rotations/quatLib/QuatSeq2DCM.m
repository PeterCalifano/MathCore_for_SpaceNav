function o_dDCM = QuatSeq2DCM(i_dQuatRot, i_bIS_VSRPplus) %#codegen
%% PROTOTYPE
% o_dDCM = QuatSeq2DCM(i_dQuatRot, i_bIS_VSRPplus) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting attitude quaternion to DCM. Convention flag set to 1
% (JPL convention) by default. 
% Multiple quaternions are automatically converted to a DCM sequence.
% REFERENCE:
% 1) Section 2.9.3, page 45, Fundamentals of Spacecraft Attitude 
% Determination and Control, Markley, Crassidis, 2014
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dQuatRot:     [4, Nseq] Input quaternion to convert to DCM. Nseq is
%                              the number of quaternions in the sequence.
% i_bIS_VSRPplus: [1] Boolean flag indicating the convention of the
%                     quaternion. 1: VSRPplus, 0: Hamilton.
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dDCM: [3, 3, Nseq] Output rotation matrix sequence
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-12-2023    Pietro Califano     Coded from reference.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Output array initialization
Nquats = size(i_dQuatRot, 1);
o_dDCM = coder.nullcopy(zeros(3, 3, Nquats));

for idV = 1:Nquats
    o_dDCM(:, :, idV) = Quat2DCM(i_dQuatRot(idV, :), i_bIS_VSRPplus);
end


end
