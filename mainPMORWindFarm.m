clear; clc;
close all

%% Load data
windFarmHyperParameters;

%% Wind Speed Range
LOW = 6;
HIGH = 10;

%% Sample p Values
% Grid sampling
K = 17; % Number of sampling vallues of p, the wind speed
ps = linspace(LOW, HIGH, K);

%% Approximate Parametric Moment using Input-Output Data
% Load input-output data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input-output data should have been collected by running 
% "windFarmFOMROM.slx" (with uncommented "Sunsystem_200") to collect input 
% data (Wt1) and FOM output response data (Yt) for wind speed 
% p = 6.4, 7.5, and 9.3, respectively. The wind speed can be assigned 
% through "windspeed.string1" in "windFarmHyperParameters.m". Then the
% commented lines 29-42 can be used to sampling data at user-defined time
% window.
%
% We directly load input-output data sampled at the time window and 
% parameter values described in the paper for ease to demonstration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Mik.mat % y(t, p)
load Uik.mat % Tuple (w(t), p))

% h = 40 * nu; % Width of time window
% sampling = 21000; % Interval between two sampling time instant
% offset = 0 * sampling; 
% 
% Mik = zeros(h*K, 1);
% Uik = zeros(h*K, nu+1); 
% 
% for i_p = 1:K
%     p = ps(i_p);
%     wt = wt_20WT_78910WS{i_p}; % "cell" data, include input w for K sampled values of p
%     yt = yt_20WT_78910WS{i_p}; % "cell" data, include output y(t, p) for K sampled values of p
%     Uik(h*(i_p-1)+1:h*i_p, :) = [wt(end-h*sampling+1-offset:sampling:end-offset, :), ones(h, 1) * p]; 
%     Mik(h*(i_p-1)+1:h*i_p, :) = yt(end-h*sampling+1-offset:sampling:end-offset, :); 
% end

% Data-driven basis function approach
% Parameters of basis functions
N = 40; % Number of basis functions

W_max = max(Uik); 
W_min = min(Uik);
W_between = W_max - W_min; 

centers = zeros(N, nu + 1); 
sigmas = ones(N, 1) * 1;
for i_b = 1:N
    for j_b = 1:nu+1
        centers(i_b, j_b) = rand() * 3 * W_between(j_b) + W_min(j_b) - W_between(j_b);
    end 
end

Rik = RBFbasisnD(Uik, centers, sigmas);
Theta_iN = Rik \ Mik;


%% Fig. 4: Output at Different Values of p
load Yt_6_4.mat
load Yt_7_5.mat
load Yt_9_3.mat
power_base = 1e6;
% Set wind speed
p = 6.4;
% Set sampling settings
offset2 = 0; 
tspan = 200; % end time of ROM
sampling2 = 5;

out = sim("windFarmFOMROM");
t = out.tout;
t = t(1:sampling2:end-offset2);
ksit = squeeze(out.Ksait.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p]; % Tuple (ksai(t), p)
kappa_iN = RBFbasisnD(wp, centers, sigmas) * Theta_iN;

figure;
plot(Yt_6_4.Time, Yt_6_4.Data./power_base, 'b', 'LineWidth', 1.2);
hold on
plot(t, kappa_iN./power_base, 'r--', 'LineWidth', 1.2);

% Set wind speed
p = 7.5;

out = sim("windFarmFOMROM");
t = out.tout;
t = t(1:sampling2:end-offset2);
ksit = squeeze(out.Ksait.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p]; % Tuple (ksai(t), p)
kappa_iN = RBFbasisnD(wp, centers, sigmas) * Theta_iN;

plot(Yt_7_5.Time, Yt_7_5.Data./power_base, 'b', 'LineWidth', 1.2);
hold on
plot(t, kappa_iN./power_base, 'r--', 'LineWidth', 1.2);

% Set wind speed
p = 9.3;

out = sim("windFarmFOMROM");
t = out.tout;
t = t(1:sampling2:end-offset2);
ksit = squeeze(out.Ksait.Data);
ksit = ksit';
ksit = ksit(1:sampling2:end-offset2, :);

n_t = size(ksit, 1);
wp = [ksit ones(n_t, 1)*p]; % Tuple (ksai(t), p)
kappa_iN = RBFbasisnD(wp, centers, sigmas) * Theta_iN;

plot(Yt_9_3.Time, Yt_9_3.Data./power_base, 'b', 'LineWidth', 1.2);
hold on
plot(t, kappa_iN./power_base, 'r--', 'LineWidth', 1.2);
xlim([80 200]);
ylim([-2 18]);
grid on
legend("FOM", "ROM");
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 15);
ylabel('Active Power (p.u.)', 'Interpreter', 'latex', 'FontSize', 15);

% %% Fig. 5: Error between FOM and ROM Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run this section, one should uncomment "Sunsystem_200" in 
% "windFarmFOMROM.slx" to collect FOM output response data (Yt) and ROM 
% state data (Ksait) under the same time series for wind speed 
% p = 6.4, 7.5, and 9.3, respectively. The wind speed can be assigned 
% through "windspeed.string1" in "windFarmHyperParameters.m".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set wind speed
% windspeed.string1 = 6.4;
% p = 6.4;
% % Set sampling settings
% offset2 = 0; 
% tspan = 200; % end time of ROM
% sampling2 = 50;
% 
% out = sim("windFarmFOMROM");
% t = out.tout;
% t = t(1:sampling2:end-offset2);
% ksit = squeeze(out.Ksait.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% Yt = squeeze(out.Yt.Data);
% Yt = Yt';
% Yt = Yt(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p]; % Tuple (ksai(t), p)
% kappa_iN = RBFbasisnD(wp, centers, sigmas) * Theta_iN;
% 
% figure;
% plot(t, Yt - kappa_iN, 'r:', 'LineWidth', 1.2);
% hold on
% 
% % Set wind speed
% windspeed.string1 = 7.5;
% p = 7.5;
% out = sim("windFarmFOMROM");
% t = out.tout;
% t = t(1:sampling2:end-offset2);
% ksit = squeeze(out.Ksait.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% Yt = squeeze(out.Yt.Data);
% Yt = Yt';
% Yt = Yt(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p]; % Tuple (ksai(t), p)
% kappa_iN = RBFbasisnD(wp, centers, sigmas) * Theta_iN;
% 
% plot(t, Yt - kappa_iN, 'g:', 'LineWidth', 1.2);
% hold on
% 
% 
% % Set wind speed
% windspeed.string1 = 9.3;
% p = 9.3;
% out = sim("windFarmFOMROM");
% t = out.tout;
% t = t(1:sampling2:end-offset2);
% ksit = squeeze(out.Ksait.Data);
% ksit = ksit';
% ksit = ksit(1:sampling2:end-offset2, :);
% 
% Yt = squeeze(out.Yt.Data);
% Yt = Yt';
% Yt = Yt(1:sampling2:end-offset2, :);
% 
% n_t = size(ksit, 1);
% wp = [ksit ones(n_t, 1)*p]; % Tuple (ksai(t), p)
% kappa_iN = RBFbasisnD(wp, centers, sigmas) * Theta_iN;
% 
% plot(t, Yt - kappa_iN, 'b:', 'LineWidth', 1.2);
% hold on
% xlim([80 200]);
% grid on
% 
% legend("p = 6.4", "p = 7.5", "p = 9.3");
% xlabel('$t$', 'interpreter', 'latex', 'fontsize', 15);
% ylabel('Error (p.u.)', 'Interpreter', 'latex', 'FontSize', 15);