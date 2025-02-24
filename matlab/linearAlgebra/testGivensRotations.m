% TEST SETUP
rng(0);

dVec2 = randn(2,1);
dTargetMatrix = randn(4,4);



%% testGivensRotVals
% TODO
[dCosTheta, dSinTheta] = GivensRotVals(dVec2);

dGrot = [dCosTheta, dSinTheta;
        -dSinTheta, dCosTheta];

dRotatedVec2 = dGrot' * dVec2;

assert(abs(dRotatedVec2(1)) > 0)
assertDifference(dRotatedVec2(2), 0.0, 1e-9);

%% testGivensRotateRows
% Test function to eliminate the second row in ui32RowIndices
ui32RowIndices = [2,3];

% Test function to eliminate an entry rotating 2 rows
dRotatedMatrix_nullRow = GivensEliminateRow(dTargetMatrix, ui32RowIndices);

% TODO add assert

%% testGivensRotateColumns
% Test function to eliminate the second column in ui32ColsIndices
ui32ColsIndices = [2,3];

[dRotatedMatrix_nullCol] = GivensEliminateColumn(dTargetMatrix, ui32ColsIndices);

% TODO add assert

%% testGivensEliminateRow
% TODO




%% testGivensEliminateCol
% TODO


