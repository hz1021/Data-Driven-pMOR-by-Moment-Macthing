%% Signal Generator
% Compute (S, L)
[S, L] = cal_QR(eigenS, eigenS_multiplicity);
L = L';
nu = size(S, 1); % The order of the reduced-order model