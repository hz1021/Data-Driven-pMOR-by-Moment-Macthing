%% N-th Approximate Parametric Moment via Data-Driven Basis Functions

%% Collect Time-Domain Input-Output Data - Training
sysg = ss(S, [], L, []); % Signal generator (S, L)
t = 0:1e-4:30; % Simulation time history
[u, ~, w] = initial(sysg, w0, t); % dim(w) = 6000 x 32;

%% Sample Y_i(p_k) for p_k's
K = 10; % num of sampling parameter points p_k
LOW = 0.1; HIGH = 1;
p_k = linspace(LOW, HIGH, K); % Sampleing set {p_k}_{i=1}^K

% Y_i = cell(K, 1);
% % Output of the full-order model
% for kk = 1:K
%     p = p_k(kk);
%     A = p*Ae + A0;
%     sysp = ss(A, B, C, []);
%     [y, ~] = lsim(sysp, u, t, x0);
%     Y_i{kk} = y;
% end

load("Y_i_16.mat");


%% Construct U_i and O_{i,K}
h = 2 * nu; % Width of time window
sampling = 410; % Interval between two sampling time instants
offset = 100000;

U_i = w(end-h*sampling+1-offset:sampling:end-offset, :);

O_iK = zeros(h*K, 1);
for k = 1:K
    O_iK((k-1)*h+1:1:(k*h), :) = Y_i{k}(end-h*sampling+1-offset:sampling:end-offset, :); % = Y_i{k}
end

%% Data-Driven: Compute Weight Matrix \widetilde\Gamma_{i,N}
%%%%% Polynomial approximation %%%%%
tic
Upsilon_K = polyBasis(p_k', N);
U_iK = kron(Upsilon_K, U_i);
vec_Gamma_iN = U_iK \ O_iK;
Gamma_iN = reshape(vec_Gamma_iN, [nu, N])';
toc
%%%%% RBF approximation %%%%%
% centers = linspace(LOW, HIGH, N)';
% sigmas = 1 * ones(N, 1); 
% Upsilon_K = RBFBasis(p_k', centers, sigmas);
% U_iK = kron(Upsilon_K, U_i);
% vec_Gamma_iN = U_iK \ O_iK;
% Gamma_iN = reshape(vec_Gamma_iN, [nu, N])';

%% Model-Based: Compute Weight Matrix \widetilde\Gamma_{N}
R_K = zeros(K, nu);
tic
for i = 1:1:K
    A = p_k(i) * Ae + A0;
    Pi_overline = sylvester(A, -S, -B*L);
    R_K(i, :) = C * Pi_overline;
end
Upsilon_K = polyBasis(p_k', N);

Gamma_N = Upsilon_K \ R_K;
toc

%% Obtain Approximate Parametric Moment \widetilde{CPi}_{i,N}(p) = Phi_N(p) * \widetilde\Gamma_{i,N}
Phi_N = cell(1, N);
for i = 1:N
   Phi_N{i} = @(p) p.^(i - 1);
   % Phi_N{i} = @(p) exp(-(p - centers(i)).^2./(2*sigmas(i)^2));
end
evaluate_Phi = @(p) cellfun(@(f) f(p), Phi_N); % Function to evaluate the value of Phi_N at specific p_k
CPI_iN_tilde = @(p) evaluate_Phi(p) * Gamma_N; % or * Gamma_iN for data-driven