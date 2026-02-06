%% N-th Approximate Parametric Moment via Series Expansion
A0_hat = A0 + esp_point*Ae;
A1_hat = Ae;

P_p_co = cell(M,1);
P_p_co{1} = sylvester(A0_hat', A0_hat, -eye(n));

for k = 2:1:M
    P_p_co{k} = sylvester(A0_hat', A0_hat, -(A1_hat' * P_p_co{k-1} + P_p_co{k-1} * A1_hat));
    if k == M
        break;
    end
end

% Pp = P_p{1} + p.*P_p{2} + p^2.*P_p{3} + p^3.*P_p{4} + h.o.t.
P_p_fun = @(p)(zeros(n,n));
for i = 1:1:M
    P_p_fun = @(p)(P_p_fun(p) + P_p_co{i}*(p-esp_point)^(i-1));
end
