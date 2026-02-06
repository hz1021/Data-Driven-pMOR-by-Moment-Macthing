%% Global POD for Parametric Model Reduction
clear; clc;

%% Full-Order Model
systemMatrices;
n = size(A0, 1); % Order of the full-order model

%% Construct the Signal Generator (S, L)
interPointsDDBF; % Interpolation points
signalGenerator; % Construct matrices S and L

% %% Collect Time-Domain Input-Output Data - Training
% sysg = ss(S, [], L, []); % Signal generator (S, L)
% t = 0:1e-4:30; % Simulation time history
% [u, ~, w] = initial(sysg, w0, t); % dim(w) = 6000 x 32;
% 
% %% Collect Snapshots X_i(p_k) for p_k's
% K = 10; % num of sampling parameter points p_k
% LOW = 0.1; HIGH = 1;
% p_k = linspace(LOW, HIGH, K); % Sampleing set {p_k}_{i=1}^K
% 
% X_i = cell(K, 1);
% % State trajectory of the full-order model
% for kk = 1:K
%     p = p_k(kk);
%     A = p*Ae + A0;
%     sysp = ss(A, B, C, []);
%     [y, ~, x] = lsim(sysp, u, t, x0);
%     X_i{kk} = x;
% end
% 
% %% Construct Global Snapshots X = [x(t_0, p_1), ..., x(t_J, p_k)] of dim n x h (~= n x K(J+1))
% h = 2 * nu; % Width of time window
% sampling = 410; % Interval between two sampling time instants
% offset = 100000;
% 
% X = zeros(n, 1);
% for k = 1:K
%     X(:, (k-1)*h+1:1:(k*h)) = X_i{k}(end-h*sampling+1-offset:sampling:end-offset, :)';
% end


%% Compute POD Modes via SVD
load("X.mat");

% POD_r(X) is by taking the r leftmost columns of the left singular matrix
tic
nu = 32; % Order of the reduced-order model, decided by the users
[V1, S, ~] = svd(X);
V = V1(:, 1:nu);
W = V;
toc


%% Construct Reduced-order Model
sysp = @(p)ss(p.*Ae + A0, B, C, 0);
sysprPOD = @(p)ss(W'*(p.*Ae + A0)*V, W'*B, C*V, 0);

%%%%% POD %%%%%
%%
p0 = 0.1;
p1 = 0.9;

xlim_min = 0; xlim_max = 3; % B. xlim_min = -2; xlim_max = 1;
xlim_bode = logspace(xlim_min, xlim_max, 1000);

[mag_FOM, ~, wout_FOM] = bode(sysp(p0), xlim_bode);
[mag_pROM, ~, wout_pROM] = bode(sysprPOD(p0), xlim_bode);

figure();
semilogx(wout_FOM, squeeze(mag_FOM),'LineWidth', 1.6, 'Color', 'b', 'LineStyle', '-');
hold on;
semilogx(wout_pROM, squeeze(mag_pROM),'LineWidth', 1.6, 'Color', 'r', 'LineStyle', '--');

[mag_FOM,phase_FOM, wout_FOM] = bode(sysp(p1), xlim_bode);
[mag_pROM,phase_pROM, wout_pROM] = bode(sysprPOD(p1), xlim_bode);

semilogx(wout_FOM, squeeze(mag_FOM),'LineWidth', 1.6, 'Color', 'b', 'LineStyle', '-');
semilogx(wout_pROM, squeeze(mag_pROM),'LineWidth', 1.6, 'Color', 'r', 'LineStyle', '--');
hold off;

title('Bode Diagram', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Magnitude (dB)', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Frequency $\omega$', 'Interpreter', 'latex', 'FontSize', 18);
legend('Original System', 'pROM with $G(p)$', 'interpreter', 'latex');
xlim([1e0 1e3]);
ylim([0.5 2.8]);
grid on;

%% Check Stability
interval = 10;
stabilityPROM = zeros(interval, 1);
for i = 1:1:interval
    p0 = 0.1 + (i-1)./interval;
    stabilityPROM(i) = isstable(sysprPOD(p0));
end
if all(stabilityPROM)
    disp('All ROMs stable')
else
    disp('Some unstable')
end