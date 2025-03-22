function dQuatSeq = qChangeConv(dQuatSeq, bH2JPL) %#codegen
%% PROTOTYPE
% [dQuatSeq] = qChangeConv(dQuatSeq, bH2JPL)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Convention switching function for quaternion array sequences. Default
% behaviour: no bH2JPL specified is equivalent to bH2JPL == true.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dQuatSeq: [4xN]   Array of N quaternions in input convention
% bH2JPL    [bool]  Flag determining the convention direction
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dQuatSeq: [4xN] Array of N quaternions in output convention
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
    bH2JPL = true;
end

% Determine array sizes
[nRows, nEntries] = size(dQuatSeq);

if nRows == 4
    % Copy array
    dQuatSeq = dQuatSeq;
elseif nEntries == 4
    % Transpose array
    dQuatSeq = dQuatSeq';
    nEntries = size(dQuatSeq, 2);
else
    assert(nRows == 4 || nEntries == 4, 'Quaternions array size is NOT consistent.')
end

% Initialize output array
dQuatSeq = zeros(4, nEntries);

if bH2JPL == true
    % Move scalar first to scalar last
    dQuatSeq(1:3, :) = dQuatSeq(2:4, :);
    dQuatSeq(4, :) = dQuatSeq(1, :);

elseif bH2JPL == false
    % Move scalar last to scalar first
    dQuatSeq(2:4, :) = dQuatSeq(1:3, :);
    dQuatSeq(1, :) = dQuatSeq(4, :);
else

end
end

