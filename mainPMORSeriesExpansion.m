clear;
clc;

%% Full-Order Model
systemMatrices;

n = size(A0, 1);

%% Construct the Signal Generator (S, L)
interPointsSE;
% Signal generator, i.e., S and L
signalGenerator;

%% Construct Reduced-Order Model
N = 4; % num of terms in the expanded series
MM = 7; % num of basis functions
M = N; % M >= N, num of terms in G(p)
esp_point = 0.55; % expansion point
% Make ROM parametric, give value of 'p = e - esp_point'
% p = e - esp_point; % Series Expansion
tic
appParaMSE;
% appParaMDDBF;
appParaPSE;
% Regularization
epsilon = 1e-14; % 1e-14
G = @(p)(PI_p_fun(p)' * P_p_fun(p) * PI_p_fun(p) + epsilon .* eye(nu))\(PI_p_fun(p)' * P_p_fun(p) * B);
toc

% G_fix = G(esp_point);
load("G_fix.mat");

sysp = @(p)ss(p.*Ae + A0, B, C, 0);

syspr = @(p)ss(S - G(p)*L, G(p), CPI_N_hat(p), 0);

sysprPOD = @(p)ss(S - G(p)*L, G(p), CPI_N_hat(p), 0);

% sysprDD = @(p)ss(S-G(p)*L, G(p), CPI_iM_tilde(p), 0);
% sysprDD = @(p)ss(S-G_fix*L, G_fix, CPI_iM_tilde(p), 0);

% sysprFix = @(p)ss(S - G_fix*L, G_fix, CPI_N_hat(p), 0);

%%%%%% Series Expansion %%%%%
%%
p0 = 0.1;
p1 = 0.9;

xlim_min = 0; xlim_max = 3; % B. xlim_min = -2; xlim_max = 1;
xlim_bode = logspace(xlim_min, xlim_max, 1000);

[mag_FOM, ~, wout_FOM] = bode(sysp(p0), xlim_bode);
[mag_pROM, ~, wout_pROM] = bode(syspr(p0), xlim_bode);

figure();
semilogx(wout_FOM, squeeze(mag_FOM),'LineWidth', 1.6, 'Color', 'b', 'LineStyle', '-');
hold on;
semilogx(wout_pROM, squeeze(mag_pROM),'LineWidth', 1.6, 'Color', 'r', 'LineStyle', '--');
% semilogx(wout_pfROM, squeeze(mag_pfROM),'LineWidth', 1.6, 'Color', '#EDB120', 'LineStyle', ':');

[mag_FOM,phase_FOM, wout_FOM] = bode(sysp(p1), xlim_bode);
[mag_pROM,phase_pROM, wout_pROM] = bode(syspr(p1), xlim_bode);

semilogx(wout_FOM, squeeze(mag_FOM),'LineWidth', 1.6, 'Color', 'b', 'LineStyle', '-');
semilogx(wout_pROM, squeeze(mag_pROM),'LineWidth', 1.6, 'Color', 'r', 'LineStyle', '--');
hold off;

stabilityFOM = isstable(sysp(p0));
% stabilityPara0 = isstable(syspr(p0));
% stabilityPara1 = isstable(syspr(p1));

title('Bode Diagram', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Magnitude (dB)', 'Interpreter', 'latex', 'FontSize', 18);
xlabel('Frequency $\omega$', 'Interpreter', 'latex', 'FontSize', 18);
% legend('Original System', 'pROM with $G(p)$', 'pROM with $\bar{G}$', 'interpreter', 'latex');
legend('Original System', 'pROM with $G(p)$', 'interpreter', 'latex');
xlim([1e0 1e3]); % B. xlim([1e-2 1e1])
ylim([0.5 2.8]); % B. ylim([-40 0]);
% ylim([0 0.6]); % for N = 100
grid on;
grid on;


%%%%% Data-Driven Basis Functions %%%%%
% %%
% p0 = 0.1;
% p1 = 0.9;
% 
% xlim_min = 0; xlim_max = 3; % B. xlim_min = -2; xlim_max = 1;
% xlim_bode = logspace(xlim_min, xlim_max, 1000);
% 
% [mag_FOM, ~, wout_FOM] = bode(sysp(p0), xlim_bode);
% [mag_pROM, ~, wout_pROM] = bode(sysprDD(p0), xlim_bode);
% 
% figure();
% semilogx(wout_FOM, squeeze(mag_FOM),'LineWidth', 1.6, 'Color', 'b', 'LineStyle', '-');
% hold on;
% semilogx(wout_pROM, squeeze(mag_pROM),'LineWidth', 1.6, 'Color', 'r', 'LineStyle', '--');
% % semilogx(wout_pfROM, squeeze(mag_pfROM),'LineWidth', 1.6, 'Color', '#EDB120', 'LineStyle', ':');
% 
% [mag_FOM,phase_FOM, wout_FOM] = bode(sysp(p1), xlim_bode);
% [mag_pROM,phase_pROM, wout_pROM] = bode(sysprDD(p1), xlim_bode);
% 
% semilogx(wout_FOM, squeeze(mag_FOM),'LineWidth', 1.6, 'Color', 'b', 'LineStyle', '-');
% semilogx(wout_pROM, squeeze(mag_pROM),'LineWidth', 1.6, 'Color', 'r', 'LineStyle', '--');
% hold off;
% 
% stabilityFOM = isstable(sysp(p0));
% % stabilityPara0 = isstable(syspr(p0));
% % stabilityPara1 = isstable(syspr(p1));
% 
% title('Bode Diagram', 'Interpreter', 'latex', 'FontSize', 18);
% ylabel('Magnitude (dB)', 'Interpreter', 'latex', 'FontSize', 18);
% xlabel('Frequency $\omega$', 'Interpreter', 'latex', 'FontSize', 18);
% % legend('Original System', 'pROM with $G(p)$', 'pROM with $\bar{G}$', 'interpreter', 'latex');
% legend('Original System', 'pROM with $G(p)$', 'interpreter', 'latex');
% xlim([1e0 1e3]); % B. xlim([1e-2 1e1])
% ylim([0.5 2.8]); % B. ylim([-40 0]);
% % ylim([0 0.6]); % for N = 100
% grid on;
% grid on;
%% Check Stability
interval = 10;
stabilityPara0v = zeros(interval, 1);
for i = 1:1:interval
    p0 = 0.1 + (i-1)./interval;
    stabilityPara0v(i) = isstable(syspr(p0));
end
if all(stabilityPara0v)
    disp('All ROMs stable')
else
    disp('Some unstable')
end
