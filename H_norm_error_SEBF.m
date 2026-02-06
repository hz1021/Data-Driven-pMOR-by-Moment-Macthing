%% Calculate Relative H2 Norm of the Error System on the Parameter Space

% Load full-order model, signal generator, and hyper parameters
mainPMORSeriesExpansion;
% mainPMORBasisFunction;
% mainPMORPOD;
% mainPMORMultiP;
% mainPMORLoewner;
% mainPMORStableMatrixInterpolation;


% Discretization of the parameter space, [p1, p2, ..., pK_bar]
K_bar = 400;
LOW = 0.1; HIGH = 1;
psn = linspace(LOW, HIGH, K_bar);

% Initialise error vectors
H2_rel_p_POD = zeros(K_bar, 1);

syserr_POD = @(p)(sysp(p) - sysprPOD(p));

for ie = 1:K_bar
    ee = psn(ie);
    % Evaluate the error system
    H2_rel_p_POD(ie) = norm(syserr_POD(ee))/norm(sysp(ee)); % Relative H2 norm of the error system
end


%% Plot the Figure
figure();
load("H2NormDD.mat"); % !!! H2NormDD should be deleted from the folder finally !!!
plot(psn, H2_rel_p_POD, linspace(LOW, HIGH, 400), H2NormDD, 'LineWidth', 1.6);
% plot(psn, H2_rel_p_POD, 'LineWidth', 1.6);
legend("POD", 'fontsize', 13.5);
set(gca, 'YScale', 'log');
grid on;
xlabel("$p$", "Interpreter", "latex", 'fontsize', 13.5);
ylabel('Relative $\mathcal{H}_2$-norm error', 'Interpreter', 'latex', 'fontsize', 13.5);