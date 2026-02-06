%% N = 1000
load("A0.mat"); load("Ae.mat"); load("B.mat"); load("C.mat");

% %% N = 100
% load("A0_100.mat"); load("Ae_100.mat"); load("B_100.mat"); load("C_100.mat");

%% Convert the sparse matrices
A0 = full(A0); Ae = full(Ae); B = full(B); C = full(C);
% A = e*Ae + A0; % e \in [0.1 1]
