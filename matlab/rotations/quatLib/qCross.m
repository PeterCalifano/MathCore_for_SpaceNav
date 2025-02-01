function q1xq2 = qCross(q1, q2) %#codegen
arguments
    q1 (4,1)
    q2 (4,1)
end
%% PROTOTYPE
% q1xq2 = qCross(q1, q2)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Computes the quaternion product between q1 and q2 using matrix Psi(q1)
% operation: [q1 x] = [Psi(q) q]. Quaternion convention: JPL (qv, qs) Left-Handed
% REFERENCE:
% 1) Fundamentals of Spacecraft Attitude Determination and Control, Markley
% Crassidis, 2014. Section 2, page 53.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% q1: [4x1] Quaternion 1 (left side of the quat. product)
% q2: [4x1] Quaternion 2 (right side of the quat. product)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% q1xq2: [4x1] Quaternion cross product between q1 and q2
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-09-2023    Pietro Califano   Coded and validated for std quat.lib.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------

q1xq2 = [q1(4) q1(3) -q1(2) q1(1); ...
         -q1(3) q1(4) q1(1) q1(2); ...
          q1(2) -q1(1) q1(4) q1(3); ...
         -q1(1) -q1(2) -q1(3) q1(4)] * q2;


end