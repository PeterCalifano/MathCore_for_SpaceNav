function [dTargetMatrix] = GivensEliminateColumn(dTargetMatrix, ui32ColumnsIndices) %#codegen
%% SIGNATURE
% [dTargetMatrix] = GivensEliminateColumn(dTargetMatrix, ui32ColumnsIndices) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% Name1                     []
% Name2                     []
% Name3                     []
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

ui32NumOfRows = size(dTargetMatrix, 2);

if coder.target("MATLAB") || coder.target("MEX") 
    % Assert sizes against indices
    % TODO
end

for ui32IdRow = 1:ui32NumOfRows

    % Get (a,b) entries
    dVal1 = dTargetMatrix(ui32IdRow, ui32ColumnsIndices(1));
    dVal2 = dTargetMatrix(ui32IdRow, ui32ColumnsIndices(2));

    % Compute Givens rotation values
    [dCos, dSin] = GivensRotVals([dVal1; dVal2]);
    
    % Rotate dTargetMatrix, equivalent to Arot = A * G
    dTargetMatrix(ui32IdRow, ui32ColumnsIndices(1)) = dCos * dVal1 - dSin * dVal2;
    dTargetMatrix(ui32IdRow, ui32ColumnsIndices(2)) = dSin * dVal1 + dCos * dVal2;

end


end
