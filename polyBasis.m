function [vector] = polyBasis(x, order)
    n = size(x, 1); % e.g., [x1; x2; x3]
    vector = zeros(n, order);
    for i = 1:order
        vector(:, i) = x.^(i-1);
    end
% vector: [1 x1 x1^2 ...; 1 x2 x2^2 ...; 1 x3 x3^2 ...]
end

