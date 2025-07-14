close all;
clear; 
clc;

% TEST SETUP
dTimeStep                = 0.01;      % s
dTimeConstFirstOrderGM   = 1.0;       % s
dSigmaFirstOrderGM       = 0.1;       % rad/s

rng(1234);  % Fix seed for reproducibility

% Generate dummy trajectory in quaternion space
dQuat0 = [1;0;0;0];
dQuat1 = [0;0.75;1.2;0.9];
dQuat1 = dQuat1/norm(dQuat1);

dTimegrid = 0:1:1000;
dNormTimegrid_ = (dTimegrid - dTimegrid(1)) / (dTimegrid(end) - dTimegrid(1)); 

dQuatSequenceOrig = InterpolateSlerp(dQuat0, dQuat1, dNormTimegrid_);

% Perturb trajectory
dQuatSequenceOut = PerturbQuaternionSequence(dQuatSequenceOrig, ...
                                             dTimeStep, ...
                                             dTimeConstFirstOrderGM, ...
                                             dSigmaFirstOrderGM);

assert(isequal(size(dQuatSequenceOut), size(dQuatSequenceOrig)), 'Output size mismatch.');

dNormOut = vecnorm(dQuatSequenceOut, 2, 1);
assert(all(abs(dNormOut - 1) < 1e-12), 'Some output quaternions are not unit norm.');

if dSigmaFirstOrderGM > 0
    % At least one sample should differ
    assert(any(max(abs(dQuatSequenceOut - dQuatSequenceOrig),[],1) > 0), ...
        'No perturbation detected (all outputs equal inputs).');
end

% Plot for visual check
dEulerOrig = quat2eul(dQuatSequenceOrig', 'ZYX');   % [yaw pitch roll]
dEulerOut  = quat2eul(dQuatSequenceOut',  'ZYX');

dTimeVec = (0:size(dQuatSequenceOrig,2)-1) * dTimeStep;

figure;
subplot(3,1,1);
plot(dTimeVec, dEulerOrig(:,1), '-', dTimeVec, dEulerOut(:,1), '--', 'LineWidth', 1.2);
ylabel('Yaw (rad)');
legend('Orig','Perturbed');
grid on;

subplot(3,1,2);
plot(dTimeVec, dEulerOrig(:,2), '-', dTimeVec, dEulerOut(:,2), '--', 'LineWidth', 1.2);
ylabel('Pitch (rad)');
grid on;

subplot(3,1,3);
plot(dTimeVec, dEulerOrig(:,3), '-', dTimeVec, dEulerOut(:,3), '--' , 'LineWidth', 1.2);
ylabel('Roll (rad)');
xlabel('Time (s)');
grid on;

sgtitle('Original vs. Perturbed Attitude (ZYX Euler)');

