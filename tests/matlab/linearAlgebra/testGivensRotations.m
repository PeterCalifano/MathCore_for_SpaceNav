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

%% testGivensEliminateRow
% Test function to eliminate the second row in ui32RowIndices
ui32TargetSubscript = [2,3];
ui32AuxRowId = 1;

[dRotatedMatrix_nullEntry23] = GivensEliminateRow(dTargetMatrix, ...
                                     ui32TargetSubscript, ...
                                     ui32AuxRowId); 

assert( abs(dRotatedMatrix_nullEntry23(ui32TargetSubscript(1), ui32TargetSubscript(2))) <= 1.5* eps)

assert( all(dRotatedMatrix_nullEntry23(ui32TargetSubscript(1), :) ~= dTargetMatrix(ui32TargetSubscript(1), :), 'all') )
assert( all(dRotatedMatrix_nullEntry23(ui32TargetSubscript(1)-1, :) ~= dTargetMatrix(ui32TargetSubscript(1)-1, :), 'all') )

%% testGivensEliminateColumn
% Test function to eliminate ui32TargetSubscript entries rotating columns
ui32TargetSubscript = [2,3];
ui32AuxColId = 2;

[dRotatedMatrix_nullEntry23_colWise] = GivensEliminateColumn(dTargetMatrix, ...
                                                              ui32TargetSubscript, ...
                                                              ui32AuxColId); 

assert( abs(dRotatedMatrix_nullEntry23_colWise(ui32TargetSubscript(1), ui32TargetSubscript(2))) <= 1.5* eps)

assert( all(dRotatedMatrix_nullEntry23_colWise(:, ui32TargetSubscript(2)) ~= dTargetMatrix(:, ui32TargetSubscript(2)), 'all') )
assert( all(dRotatedMatrix_nullEntry23_colWise(:, ui32TargetSubscript(2)-1) ~= dTargetMatrix(:, ui32TargetSubscript(2)-1), 'all') )

%% testGivensEliminateQR
% Test function to perform QR decomposition using Givens rotations
[dR] = GivensEliminateQR(dTargetMatrix, true);

% Assert eliminated matrix is upper triangular
assert( istriu(dR) );

% TODO, requires dQ computation
% Ensure Q is orthogonal: Q'Q should be identity
% assert( norm(dQ' * dQ - eye(size(dQ, 1))) < 1e-10 );

% Ensure Q * R reconstructs the original matrix
% assert( norm(dQ * dR - dTargetMatrix) < 1e-10 );

% Check dimensions of Q and R
% [ui32NumOfRows, ui32NumOfCols] = size(dTargetMatrix);
% assert( isequal(size(dQ), [ui32NumOfRows, ui32NumOfCols]) );  % Q should be square (m x m)
% assert( isequal(size(dR), [ui32NumOfRows, ui32NumOfCols]) );  % R should have the same size as the original matrix

