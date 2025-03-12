function dScatteredValues = uniformScatter(dMinVec, dMaxVec, ui32NumPoints) %#codegen
arguments
    dMinVec         (:,1) double {isnumeric}
    dMaxVec         (:,1) double {isnumeric}
    ui32NumPoints   (1,1) uint32 {isscalar, isnumeric}
end
%% PROTOTYPE
% dScatteredValues = uniformScatter(dMinVec, dMaxVec, ui32NumPoints) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function randomly generating uniformly distributed 1D values given lower and upper bounds of the interval
% and the desired number of points.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dMinVec         (:,1) double {isnumeric}
% dMaxVec         (:,1) double {isnumeric}
% ui32NumPoints   (1,1) uint32 {isscalar, isnumeric}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dScatteredValues
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-03-2024        Pietro Califano     First version coded (scalar values)
% 11-03-2025        Pietro Califano     Update of function to new coding standards     
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Get size of inputs
ui32VecSpaceDim = size(dMinVec, 1);
assert(ui32VecSpaceDim == size(dMaxVec, 1), "Max and Min vectors have different dimensions")

dRandPoints = transpose(rand(ui32NumPoints, size(dMaxVec, 1), size(dMaxVec, 2)));
dScatteredValues = dMinVec + (dRandPoints .* (dMaxVec - dMinVec));

end

