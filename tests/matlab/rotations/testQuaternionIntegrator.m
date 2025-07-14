close all
clear
clc

% TEST SETUP
objQuatInt = QuaternionIntegrator();


%% test_LieGroupMode
fOmega = @(dT) [0; 0; deg2rad(90)]; % constant angular velocity
dQuat0 = [1; 0; 0; 0];
dQuatEndLie   = objQuatInt.Integrate(dQuat0, fOmega, 0, 2, 0.01, 'lie');
dQuatEndRKMK4 = objQuatInt.Integrate(dQuat0, fOmega, 0, 2, 0.01, 'rkmk4');

disp('Lie-Euler result:');
disp(dQuatEndLie);
disp('RKMK4 result:');
disp(dQuatEndRKMK4);

%% test_RKMK4Mode
