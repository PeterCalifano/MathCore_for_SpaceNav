function o_dQuatRot = DCM2quatSeq(i_dDCM, i_bIS_VSRPplus) %#codegen
%% PROTOTYPE
% o_dQuatRot = DCM2quat(i_dDCM, i_bIS_VSRPplus)%#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function converting a DCM to Attitude quaternion. Conversion occurs
% according to JPL convention. Set i_bIS_JPL_CONV = false if Hamilton
% convention is being used. The conversion is made numerically optimal.
% Multiple DCM are automatically converted to a quaternion sequence.
% REFERENCE:
% 1) Section 2.9.3, page 48, Fundamentals of Spacecraft Attitude 
% Determination and Control, Markley, Crassidis, 2014
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% o_dDCM:         [3, 3, Nseq]  Input rotation matrix to convert. Nseq is
%                              the number of quaternions in the sequence. 
% i_bIS_VSRPplus: [1]  Boolean flag indicating the convention of the
%                      quaternion. 1: VSRPplus, 0: SVRPplus (as MATLAB)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dQuatRot:     [4, Nseq]  Output quaternion sequence
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 16-12-2023    Pietro Califano     Coded from reference. Validated against
%                                   MATLAB library.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Vectorization for multiple quaternions in one call (improvement using
% vectorization over third dimension).
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Quaternion output initialization
Nquats = size(i_dDCM, 3);
o_dQuatRot = coder.nullcopy(zeros(4, Nquats));

for idV = 1:Nquats
    % Compute trace of DCM
    trDCM = trace(i_dDCM(:, :, idV));

    [~, idMax] = max([i_dDCM(1,1, idV), i_dDCM(2,2, idV), i_dDCM(3,3, idV), trDCM]);

    % Adjust computation case for numerical accuracy to convert DCM to NON
    % UNITARY QUATERNION
    switch idMax
        case 4
            % Scalar part is max
            o_dQuatRot(:, idV) = [i_dDCM(2, 3, idV) - i_dDCM(3, 2, idV);
                i_dDCM(3, 1, idV) - i_dDCM(1, 3, idV);
                i_dDCM(1, 2, idV) - i_dDCM(2, 1, idV);
                1 + trDCM];
        case 3
            % Third vector component is max
            o_dQuatRot(:, idV) = [i_dDCM(3, 1, idV) + i_dDCM(1, 3, idV);
                                  i_dDCM(3, 2, idV) + i_dDCM(2, 3, idV);
                                  1 + 2*i_dDCM(3, 3, idV) - trDCM;
                                  i_dDCM(1, 2, idV) - i_dDCM(2, 1, idV)];
                  
        case 2
            % Second vector component is max
            o_dQuatRot(:, idV) = [i_dDCM(2, 1, idV) + i_dDCM(1, 2, idV);
                1 + 2*i_dDCM(2, 2, idV) - trDCM;
                i_dDCM(2, 3, idV) + i_dDCM(3, 2, idV);
                i_dDCM(3, 1, idV) - i_dDCM(1, 3, idV)];

        case 1
            % First vector component is max
            o_dQuatRot(:, idV) = [1 + 2*i_dDCM(1, 1, idV) - trDCM;
                i_dDCM(1, 2, idV) + i_dDCM(2, 1, idV);
                i_dDCM(1, 3, idV) + i_dDCM(3, 1, idV);
                i_dDCM(2, 3, idV) - i_dDCM(3, 2, idV)];
    end

end

% Normalize quaternion
o_dQuatRot = o_dQuatRot./vecnorm(o_dQuatRot, 2, 1);

% Assign quaternion components depending on quaternion convention
if i_bIS_VSRPplus == false
    % SVRPplus convention (false) --> SWAP
    o_dQuatRot = [o_dQuatRot(4, :);
                  o_dQuatRot(1, :);
                  o_dQuatRot(2, :);
                  o_dQuatRot(3, :)];

elseif i_bIS_VSRPplus == true
    % VSRPplus convention (false). DO NOTHING: DEFAULT
    %     o_dQuatRot(1) = qv1;
    %     o_dQuatRot(2) = qv2;
    %     o_dQuatRot(3) = qv3;
    %     o_dQuatRot(4) = qs ;
end

end


