%% Multi Para MM for Parametric Model Reduction
clear; clc;

%% Full-Order Model
systemMatrices;
n = size(A0, 1); % Order of the full-order model

%% Construct the Signal Generator (S, L)
% interPointsDDBF; % Interpolation points
interPointsSE; % Interpolation points
signalGenerator; % Construct matrices S and L

%% Compute Projection Matrix via Arnoldi
tic
np = 0; % The matched order of the parameter p
shiftPoles = eigenS;
% shiftPoles = 0;
V = [];
for si = 1:1:length(eigenS)
    [Vi, ~] = arnoldi((shiftPoles(si)*eye(n)-A0-0.55*Ae)\Ae, (shiftPoles(si)*eye(n)-A0-0.55*Ae)\B, np);
    V = [V, real(Vi), imag(Vi)];
end
[V, ~] = qr(V, 0); % Economy-size QR
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