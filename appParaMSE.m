%% N-th Approximate Parametric Moment via Series Expansion
A0_hat = A0 + esp_point*Ae;
A1_hat = Ae;

PI_p_co = cell(M,1);
PI_p_co{1} = sylvester(A0_hat,-S,-B*L);

for k = 2:1:N
    PI_p_co{k} = sylvester(A0_hat,-S,-A1_hat*PI_p_co{k-1});
end
for k = N+1:1:M
    PI_p_co{k} = zeros(n, nu);
end

% PIp = PI_p{1} + p.*PI_p{2} + p^2.*PI_p{3} + p^3.*PI_p{4} + h.o.t.
PI_p_fun = @(p)(zeros(n,nu));
for i = 1:1:N
    PI_p_fun = @(p)(PI_p_fun(p) + PI_p_co{i}*(p-esp_point)^(i-1));
end

CPI_N_hat = @(p)(C*PI_p_fun(p));
