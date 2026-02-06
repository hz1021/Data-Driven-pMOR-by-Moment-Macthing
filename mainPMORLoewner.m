%% Para Loewner Framework for Parametric Model Reduction
clear; clc;

%% Full-Order Model
systemMatrices;
n = size(A0, 1); % Order of the full-order model

%% Construct the Signal Generator (S, L)
interPointsDDBF; % Interpolation points, interPointsDDBF
signalGenerator; % Construct matrices S and L

%% Collect Frequency-Domain Data - Training
% Interpolation points
lambda = eigenS;
mu = - eigenS;
% Dimension
nL = length(lambda);
pL = length(mu);

K = 10; % num of sampling parameter points p_k
LOW = 0.1; HIGH = 1;
p_k = linspace(LOW, HIGH, K); % Sampleing set {p_k}_{i=1}^K
% Sampling parameter values
pii = p_k(1:2:end);
nuL = p_k(2:2:end);
% Dimension
mL = length(pii);
rL = length(nuL);

%% Collect Data for Loewner Matrix
sysp = @(p)ss(p.*Ae + A0, B, C, 0);
w = zeros(nL, mL);
for in = 1:1:nL
    for im = 1:1:mL
        w(in, im) = freqresp(sysp(pii(im)), lambda(in));
    end
end

v = zeros(pL, rL);
for ip = 1:1:pL
    for ir = 1:1:rL
        v(ip, ir) = freqresp(sysp(nuL(ir)), mu(ip));
    end
end

tic
L = zeros(pL*rL, nL*mL);
for i = 1:1:nL
    for k = 1:1:pL
        ell = zeros(rL, mL);
        for j = 1:1:mL
            for l = 1:1:rL
                ell(l, j) = (v(k, l) - w(i, j)) / ((mu(k) - lambda(i)) * (nuL(l) - pii(j)));
            end
        end
        L(1+(k-1)*rL:k*rL, 1+(i-1)*mL:i*mL) = ell;
    end
end

%% Construct \mathbb{A} and \mathbb{B} Matrices
[~, SL, VL] = svd(L);
c = VL(:, end);
AM = reshape(c, [mL, nL]);

BM = zeros(mL, nL);
for i = 1:1:mL
    for j = 1:1:nL
        BM(i, j) = AM(i, j) * w(j, i);
    end
end

pV = zeros(mL, 1);
for i = 1:1:mL
    pV(i) = 1./ (prod(pii) ./ pii(i));
end


%% Construct Reduced-order Model
E = zeros((nL-1+mL+mL), (nL+mL-1+mL));
E(1:(nL-1), 1:nL) = [ones(nL-1, 1), -eye(nL-1)];

Arom1 = [-lambda(1)*ones(nL-1, 1), diag(lambda(2:end))];
Arom2 = @(p)[(p-pii(1))*ones(1, mL-1); diag(pii(2:end) - p)];
Arom = @(p)[Arom1, zeros(nL-1, mL-1), zeros(nL-1, mL);
    AM, Arom2(p), zeros(mL, mL);
    BM, zeros(mL, mL-1), [Arom2(p), pV]];

Brom = [zeros(nL-1, 1); pV; zeros(mL, 1)];
Crom = zeros(1, nL+mL+mL-1);
Crom(end) = -1;
toc

sysprPOD = @(p)dss(-Arom(p), Brom, Crom, 0, E);

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