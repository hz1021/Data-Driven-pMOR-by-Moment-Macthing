function [V,H] = arnoldi(A, b, np)
n = length(A);
V = zeros(n, np);
b = b/norm(b);
V(:, 1) = b;
H = zeros(min(np+1, np), n);
for k = 1:np
    w = A * V(:,k);
    for i = 1:k
        H(i, k) = V(:, i)' * w;
        w = w - H(i, k) * V(:, i);
    end
    if k < n
       H(k+1, k) = norm(w);
       if H(k+1, k) == 0, return, end
       V(:, k+1) = w / H(k+1, k);
   end
end