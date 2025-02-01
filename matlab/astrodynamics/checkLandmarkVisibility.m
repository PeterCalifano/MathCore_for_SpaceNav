function [outputArg1,outputArg2] = checkLandmarkVisibility(inputArg1,inputArg2)
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% Name4                     []
% Name5                     []
% Name6                     []
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% DD-MM-YYYY        Pietro Califano         Modifications
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

%% DEVNOTE Options:
%  1) Assume cone-shaped Field of View and check if the angle between the LoS to the landmark and the camera
%  boresight is < threshold
%  2) Consider the camera reference frame axis and compute the spherical angles of the line-of-sight,
%  compare them with the maximum angles determined by the Field of View
%  3) Project the point in the detector place/image plane and checl the pixel coordinates against the 
%  boundaries of the image. The landmark is not visible if outside, visible if inside. This option can more
%  easily account for lens distorsions.


end