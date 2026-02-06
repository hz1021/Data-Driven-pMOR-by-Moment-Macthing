%% Matrix Interpolation BT for Stability-Preserving Parametric Model Reduction
clear; clc;

%% Full-Order Model
systemMatrices;
n = size(A0, 1); % Order of the full-order model

%% Step A. Compute Reduced Matrices for Each Parameter Point
K = 10; % num of sampling parameter points p_k
LOW = 0.1; HIGH = 1;
p_k = linspace(LOW, HIGH, K); % Sampleing set {p_k}_{i=1}^K

sysp = @(p)ss(p.*Ae + A0, B, C, 0);

Ari = cell(K, 1);
Bri = cell(K, 1);
Cri = cell(K, 1);
nu = 32; % Reduced order

%% Choice 1: BT
tic
Vi = cell(K, 1);
Wi = cell(K, 1);

for i = 1:1:K
    [~, ~, TL, TR] = balreal(sysp(p_k(i)));
    TR = TR(:, 1:nu);
    TL = TL(1:nu, :);

    Vi{i} = TR(:, 1:nu);
    Wi{i} = TL(1:nu, :);

    Ari{i} = TL * sysp(p_k(i)).A * TR;
    Bri{i} = TL * sysp(p_k(i)).B;
    Cri{i} = sysp(p_k(i)).C * TR;
end

% %% Choice 2: MM
% interPointsDDBF; % Interpolation points
% signalGenerator; % Construct matrices S and L
% 
% % Compute Projection Matrix via Arnoldi
% np = 0; % The matched order of the parameter p
% shiftPoles = eigenS;
% Vi = cell(K, 1);
% 
% for i = 1:1:K
%     V = [];
%     for si = 1:1:length(eigenS)
%         [Vitemp, ~] = arnoldi((shiftPoles(si)*eye(n)-A0-p_k(i)*Ae)\eye(n), (shiftPoles(si)*eye(n)-A0-p_k(i)*Ae)\B, np);
%         V = [V, real(Vitemp), imag(Vitemp)];
%     end
%     [V, ~] = qr(V, 0); % Economy-size QR
% 
%     Vi{i} = V;
% 
%     Ari{i} = V' * sysp(p_k(i)).A * V;
%     Bri{i} = V' * sysp(p_k(i)).B;
%     Cri{i} = sysp(p_k(i)).C * V;
% end


%% Step B. Transformation for Stability
Li = cell(K, 1);
for i = 1:1:K
    Qi = Cri{i}' * Cri{i};
    Pri = lyap(Ari{i}', Ari{i}, Qi);
    Li{i} = chol(Pri);
end

%% Step C. Transformation for Generalised Coordinates
V0 = Vi{1};

VallM = [];
for i = 1:1:K
    VallM = [VallM, Vi{i}];
end

[R, ~, ~] = svd(VallM);
V0 = R(:, 1:nu);

Ri = cell(K, 1);
for i = 1:1:K
    Si = (Vi{i} / Li{i})' * V0;
    [Ui, ~, Zi] = svd(Si);
    Ri{i} = Ui * Zi';
end

%% Step D. Interpolation Process
Ari_hat = cell(K, 1);
Bri_hat = cell(K, 1);
Cri_hat = cell(K, 1);
for i = 1:1:K
    Ari_hat{i} = Ri{i}' * Li{i} * Ari{i} / Li{i} * Ri{i};
    Bri_hat{i} = Ri{i}' * Li{i} * Bri{i};
    Cri_hat{i} = Cri{i} / Li{i} * Ri{i};
end

% Construct reduced-order model
Ar = @(p)zeros(nu, nu);
Br = @(p)zeros(nu, 1);
Cr = @(p)zeros(1, nu);

deltaP = 0.1;
for i = 1:1:K
    Ar = @(p) Ar(p) + max((1 - abs(p - p_k(i)) ./ deltaP), 0) * Ari_hat{i};
    Br = @(p) Br(p) + max((1 - abs(p - p_k(i)) ./ deltaP), 0) * Bri_hat{i};
    Cr = @(p) Cr(p) + max((1 - abs(p - p_k(i)) ./ deltaP), 0) * Cri_hat{i};
    % Ar = @(p) Ar(p) + 0.1 * Ari_hat{i};
    % Br = @(p) Br(p) + 0.1 * Bri_hat{i};
    % Cr = @(p) Cr(p) + 0.1 * Cri_hat{i};
end
toc
sysprPOD = @(p) ss(Ar(p), Br(p), Cr(p), 0);



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