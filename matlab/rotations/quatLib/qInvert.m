function o_dInvQuatSeq = qInvert(i_dQuatSeq, i_bIS_VSRPplus) %#codegen
arguments
    i_dQuatSeq     (4, :)
    i_bIS_VSRPplus (1,1) logical
end
%% PROTOTYPE
% o_dInvQuatSeq = qInvert(i_dQuatSeq, i_bIS_VSRPplus)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the inverse according to VSRP+ convention (right-handed). Set i_bIS_VSRPplus = false if SVRP+ 
% convention with scalar first is being used. 
% (SV) Scalar first, Vector last
% (P) Passive 
% (R) Successive coordinate transformations have the unmodified quaternion chain on the Right side of
%     the triple product.
% (plus) Right-Handed Rule for the imaginary numbers i, j, k. (aka Hamilton)
% REFERENCE:
% 1) Fundamentals of Spacecraft Attitude Determination and Control, Markley
% Crassidis, 2014. Section 2.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dQuatSeq:     [4, N] Sequence of attitude quaternion with convention 
%                        determined by i_bConvFlag. 
% i_bIS_VSRPplus:   [1]  Boolean flag indicating the convention of the
%                       quaternion. 1: VSRPplus, 0: SVRPplus (as MATLAB)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dInvQuatSeq: [4, N] Sequence of inverted attitude quaternions
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 25-11-2023    Pietro Califano   Function coded.
% 14-12-2023    Pietro Califano   Change of convention flag for more
%                                 self-explainable code.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
if nargin < 2
    % Default behaviour
    i_bIS_VSRPplus = true;
end

assert(islogical(i_bIS_VSRPplus), 'Function expected a boolean flag as second input, but found a different dtype!')

if i_bIS_VSRPplus == true
    o_dInvQuatSeq = [-i_dQuatSeq(1:3, :); i_dQuatSeq(4, :)];
else
    o_dInvQuatSeq = [i_dQuatSeq(1, :); -i_dQuatSeq(2:4, :)];
end

end
