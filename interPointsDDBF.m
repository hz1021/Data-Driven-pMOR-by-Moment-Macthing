%% Selected Interpolation Points in Illustration of Basis Function
s = linspace(0, 3.1, 16); % 16 equaldistant sampling points on [0, 3.1]
eigenS = 10.^s.*1i; 
eigenS_multiplicity = 1.*ones(1, size(eigenS, 2));

nu = sum(eigenS_multiplicity) * 2; % Order of the reduced-order model


%% Initial Conditions for data-driven settings
% Initial value of \omega, i.e., the signal generator
w0_value = 0.3;
rng(1)
w0 = w0_value * rand(nu, 1); 

% Initial value of the full-order model
x0_value = 0;
rng(2)
x0 = x0_value * rand(n, 1); 