function scatteredValues = uniformScatter(MinVec, MaxVec, nPoints) 
%% PROTOTYPE
% scatteredValues = uniformScatter(MinVec, MaxVec, nPoints) 
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function randomly generating uniformly distributed 1D values given lower and upper bounds of the interval
% and the desired number of points.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% MinVec
% MaxVec
% nPoints
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% scatteredValues
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 05-03-2024        Pietro Califano         First version coded (scalar values)
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% 1) Version supporting input vectors for uniform scattering
% -------------------------------------------------------------------------------------------------------------
%% Function code
% Get size of inputs
vecSpaceDim = size(MinVec, 1);
assert(vecSpaceDim == size(MaxVec, 1), "Max and Min vectors have different dimensions")

if vecSpaceDim == 1
    % 1D uniform random scatter
    randPoints = rand(nPoints, 1);
    scatteredValues = MinVec + (randPoints * (MaxVec - MinVec));

else
    error('Only 1D vector spaces currently supported.')
end

end

