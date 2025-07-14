close all
clear
clc


% TEST SETUP
dTestVector = [1, 1, 0;
               0, 0, 1;
               0, 0, 0];

dDirections = [0, 0, 0;
               0, 0, 0;
               1, 1, 1];

dAngles = deg2rad([90, 0, -90]);

%% test_Rot3dVecAboutDir_fullVect
dRotateVectors = Rot3dVecAboutDir(dDirections, dTestVector, dAngles);

assertDifference(dRotateVectors(:, 1), [0;1;0], 1e-12);
assertDifference(dRotateVectors(:, 2), [1;0;0], 1e-12);
assertDifference(dRotateVectors(:, 3), [1;0;0], 1e-12);

%% test_Rot3dVecAboutDir_conveniency cases
dAngles = deg2rad(90);

dRotateVectors = Rot3dVecAboutDir(dDirections(:,1), dTestVector, dAngles);

assertDifference(dRotateVectors(:, 1), [0;1;0], 1e-12);
assertDifference(dRotateVectors(:, 2), [0;1;0], 1e-12);
assertDifference(dRotateVectors(:, 3), [-1;0;0], 1e-12);
