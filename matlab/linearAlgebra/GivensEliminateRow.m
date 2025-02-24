function [dTargetMatrix] = GivensEliminateRow(dTargetMatrix, ui32RowsIndices) %#codegen
%% SIGNATURE
% [dTargetMatrix] = GivensEliminateRow(dTargetMatrix, ui32RowsIndices) %#codegen
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
% 24-02-2025    Pietro Califano     First version implemented
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

ui32NumOfCols = size(dTargetMatrix, 2);

if coder.target("MATLAB") || coder.target("MEX") 
    % Assert sizes against indices
    % TODO
end

for ui32IdCol = 1:ui32NumOfCols

    % Get (a,b) entries
    dVal1 = dTargetMatrix(ui32RowsIndices(1), ui32IdCol);
    dVal2 = dTargetMatrix(ui32RowsIndices(2), ui32IdCol);

    % Compute Givens rotation values
    [dCos, dSin] = GivensRotVals([dVal1; dVal2]);
    
    % Rotate dTargetMatrix, equivalent to Arot = G^T * A
    dTargetMatrix(ui32RowsIndices(1), ui32IdCol) = dCos * dVal1 - dSin * dVal2;
    dTargetMatrix(ui32RowsIndices(2), ui32IdCol) = dSin * dVal1 + dCos * dVal2;

end

end











