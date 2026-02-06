%% Matrix Interpolation BT for Parametric Model Reduction
clear; clc;

%% Full-Order Model
systemMatrices;
n = size(A0, 1); % Order of the full-order model

%% Construct BT Projection Matrices
K = 10; % num of sampling parameter points p_k
LOW = 0.1; HIGH = 1;
p_k = linspace(LOW, HIGH, K); % Sampleing set {p_k}_{i=1}^K

sysp = @(p)ss(p.*Ae + A0, B, C, 0);

Vall = cell(K, 1);
Wall = cell(K, 1);
nu = 32; % Reduced order
for i = 1:1:K
    [~, ~, TL, TR] = balreal(sysp(p_k(i)));
    Vall{i} = TR(:, 1:nu);
    TLTran = TL';
    Wall{i} = TLTran(:, 1:nu);
end

for i = 1:1:K
    [~, ~, TL, TR] = balreal(sysp(p_k(i)));
    TR = TR(:, 1:nu);
    TL = TL(1:nu, :);

    Aric{i} = TL * sysp(p_k(i)).A * TR;
    Bric{i} = TL * sysp(p_k(i)).B;
    Cric{i} = sysp(p_k(i)).C * TR;
end

VallM = [];
for i = 1:1:K
    VallM = [VallM, Vall{i}];
end

[R, ~, ~] = svd(VallM);
R = R(:, 1:nu);
R = Vall{1};

%% Compute Reduced Matrices for Each Parameter Point
Ari = cell(K, 1);
Bri = cell(K, 1);
Cri = cell(K, 1);
for i = 1:1:K
    Ari{i} = (Wall{i}' * R) \ Aric{i} / (R' * Vall{i});
    Bri{i} = (Wall{i}' * R) \ Bric{i};
    Cri{i} = Cric{i} / (R' * Vall{i});
    % Ari{i} = (Wall{i}' * sysp(p_k(i)).A * Vall{i});
    % Bri{i} = (Wall{i}' * sysp(p_k(i)).B);
    % Cri{i} = (sysp(p_k(i)).C * Vall{i});
end

%% Construct Reduced-order Model
Ar = @(p)zeros(nu, nu);
Br = @(p)zeros(nu, 1);
Cr = @(p)zeros(1, nu);

for i = 1:1:K
    Ar = @(p) Ar(p) + max((1 - abs(p - p_k(i)) ./ 0.1), 0) * Ari{i};
    Br = @(p) Br(p) + max((1 - abs(p - p_k(i)) ./ 0.1), 0) * Bri{i};
    Cr = @(p) Cr(p) + max((1 - abs(p - p_k(i)) ./ 0.1), 0) * Cri{i};
end
sysprPOD = @(p) ss(Ar(p), Br(p), Cr(p), 0);

% i = 9;
% Ari{1} = (Wall{i}' * sysp(p_k(i)).A * Vall{i});
% Bri{1} = (Wall{i}' * sysp(p_k(i)).B);
% Cri{1} = (sysp(p_k(i)).C * Vall{i});
% sysprPOD = @(p) ss(Ari{1}, Bri{1}, Cri{1}, 0);

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