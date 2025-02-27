%% Some quick experiments with GPT below
% Setup
M  = 5;    % number of measurements (rows)
Nx = 3;    % dimension of camera/IMU state
Nf = 2;    % dimension of feature state

% Example: let Hf have some partial linear dependence on Hx
Hx = randn(M, Nx);
Hf = 0.5 * Hx(:, 1:Nf) + 0.000000001 * randn(M, Nf);

% Concatenate full Jacobian
H = [Hx, Hf];  % size M x (Nx + Nf)

% Residual vector
r = H * randn(5,1);

% Form augmented matrix [H | r]
H_aug = [H, r];  % size M x (Nx + Nf + 1)

% Test QR first
[Q, R] = qr(H_aug);

% Print
disp(Q)
disp(R)

H * Q(:, 4:5)

A = null(H, 1e-6)
H * A
A' * H'


