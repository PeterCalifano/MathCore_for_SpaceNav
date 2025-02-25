% TEST SETUP
rng(0);

dVec2 = randn(2,1);
dTargetMatrix = randn(4,4);



%% testGivensRotVals
[dCosTheta, dSinTheta] = GivensRotVals(dVec2);

dGrot = [dCosTheta, dSinTheta;
        -dSinTheta, dCosTheta];

dRotatedVec2 = dGrot' * dVec2;

assert(abs(dRotatedVec2(1)) > 0)
assertDifference(dRotatedVec2(2), 0.0, 1e-9);

%% testGivensRotateRows
% Test function to eliminate the second row in ui32RowIndices
ui32TargetSubscript = [2,3];
ui32AuxRowId = 1;

[dRotatedMatrix_nullEntry23] = GivensEliminateRow(dTargetMatrix, ...
                                     ui32TargetSubscript, ...
                                     ui32AuxRowId); %#codegen

assert( abs(dRotatedMatrix_nullEntry23(ui32TargetSubscript(1), ui32TargetSubscript(2))) <= 1.5* eps)

assert( all(dRotatedMatrix_nullEntry23(ui32TargetSubscript(1), :) ~= dTargetMatrix(ui32TargetSubscript(1), :), 'all') )
assert( all(dRotatedMatrix_nullEntry23(ui32TargetSubscript(1)-1, :) ~= dTargetMatrix(ui32TargetSubscript(1)-1, :), 'all') )

%% testGivensRotateColumns
% Test function to eliminate ui32TargetSubscript entries rotating columns
ui32TargetSubscript = [2,3];
ui32AuxColId = 2;

[dRotatedMatrix_nullEntry23_colWise] = GivensEliminateColumn(dTargetMatrix, ...
                                                              ui32TargetSubscript, ...
                                                              ui32AuxColId); 

assert( abs(dRotatedMatrix_nullEntry23_colWise(ui32TargetSubscript(1), ui32TargetSubscript(2))) <= 1.5* eps)

assert( all(dRotatedMatrix_nullEntry23_colWise(:, ui32TargetSubscript(2)) ~= dTargetMatrix(:, ui32TargetSubscript(2)), 'all') )
assert( all(dRotatedMatrix_nullEntry23_colWise(:, ui32TargetSubscript(2)-1) ~= dTargetMatrix(:, ui32TargetSubscript(2)-1), 'all') )

%% testGivensEliminateRow
% TODO




%% testGivensEliminateCol
% TODO






