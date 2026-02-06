function [vector] = RBFBasis(x, x_c, sigma)
%POLYBASIS Summary of this function goes here
%   Detailed explanation goes here
n = size(x, 1);  % [x1; x2; x3]
n_basis = size(x_c, 1); 
vector = zeros(n, n_basis);
for i = 1:n_basis
   vector(:, i) = exp(-(x - x_c(i)).^2./(2*sigma(i)^2));
end
% vector: [1 x1 x1^2 ...; 1 x2 x2^2 ...; 1 x3 x3^2 ...]
end
