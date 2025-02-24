function dTargetMatrix = GivensRotateRows(dTargetMatrix, ui32EntrySubscript) %#codegen
arguments
    dTargetMatrix       (:,:) double {ismatrix}
    ui32EntrySubscript  (1,2) uint32 {isvector}
end
%% SIGNATURE
% dTargetMatrix = GivensRotateRows(dTargetMatrix, ui32EntrySubscript) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Column-wise Givens elimination algorithm. This affets the entries of dTargetMatrix in two rows only, 
% leaving all others unchanged. It can be used to perform elimination of entries along columns.
% REFERENCE: Golub, 2013, Matrix Computations, pag 241.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dTargetMatrix       (:,:) double {ismatrix}
% ui32EntrySubscript  (1,2) uint32 {isvector}
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dTargetMatrix       (:,:) double {ismatrix}
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 24-02-2025    Pietro Califano     First version implemented.
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

% for ui32IdCol = 1:ui32NumOfRows

    % Get (a,b) entries
    dVal1 = dTargetMatrix(ui32EntrySubscript(1), ui32IdCol);
    dVal2 = dTargetMatrix(ui32EntrySubscript(2), ui32IdCol);

    % Compute Givens rotation values
    [dCos, dSin] = GivensRotVals([dVal1; dVal2]);
    
    % Rotate dTargetMatrix, equivalent to Arot = G^T * A
    dTargetMatrix(ui32EntrySubscript(1), ui32IdCol) = dCos * dVal1 - dSin * dVal2;
    dTargetMatrix(ui32EntrySubscript(2), ui32IdCol) = dSin * dVal1 + dCos * dVal2;

% end

end
