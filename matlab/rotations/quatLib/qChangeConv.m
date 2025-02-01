function o_dQuatSeq = qChangeConv(i_dQuatSeq, i_bH2JPL) %#codegen
%% PROTOTYPE
% [o_dQuatSeq] = qChangeConv(i_dQuatSeq, i_bH2JPL)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Convention switching function for quaternion array sequences. Default
% behaviour: no i_bH2JPL specified is equivalent to i_bH2JPL == true.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dQuatSeq: [4xN]   Array of N quaternions in input convention
% i_bH2JPL    [bool]  Flag determining the convention direction
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_dQuatSeq: [4xN] Array of N quaternions in output convention
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 25-11-2023    Pietro Califano     Coded for std quat. library function.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Default behaviour Hamilton to JPL conversion
if nargin < 2
    i_bH2JPL = true;
end

% Determine array sizes
[nRows, nEntries] = size(i_dQuatSeq);

if nRows == 4
    % Copy array
    dQuatSeq = i_dQuatSeq;
elseif nEntries == 4
    % Transpose array
    dQuatSeq = i_dQuatSeq';
    nEntries = size(dQuatSeq, 2);
else
    assert(nRows == 4 || nEntries == 4, 'Quaternions array size is NOT consistent.')
end

% Initialize output array
o_dQuatSeq = zeros(4, nEntries);

if i_bH2JPL == true
    % Move scalar first to scalar last
    o_dQuatSeq(1:3, :) = dQuatSeq(2:4, :);
    o_dQuatSeq(4, :) = dQuatSeq(1, :);

elseif i_bH2JPL == false
    % Move scalar last to scalar first
    o_dQuatSeq(2:4, :) = dQuatSeq(1:3, :);
    o_dQuatSeq(1, :) = dQuatSeq(4, :);
else

end
end

