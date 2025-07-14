close all
clear
clc

%% test_UniformlySampleSO3HaarDistr
N = 1000; % Number of random rotations
points = zeros(N,3);

dR = UniformlySampleSO3HaarDistr(N);

for i = 1:N
    points(i,:) = dR(:,:,i) * [1; 0; 0]; % Apply rotation to x-axis unit vector
end

% Scatter plot of rotated points
scatter3(points(:,1), points(:,2), points(:,3), 'filled');
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Random Rotations Using Haar Measure');
