%% Selected Interpolation Points in Illustration of Series Expansion
s = linspace(0, 3.1, 50); % 50 equaldistant sampling points on [0, 3.1]
eigenS = 10.^s.*1i;
eigenS_multiplicity = 1.*ones(1, size(eigenS, 2));

nu = sum(eigenS_multiplicity) * 2; % Order of the reduced-order model