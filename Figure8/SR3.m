function [Xi_full,Xi_sparse,trimmed_points] = SR3(A, dxdt,tau, kappa, ss, tf, maxiter)
% A: function library
% dxdt: derivative 
% N: number of equations
% tau: threshold for coefficients
% kappa: controls how closely sparse representation must match least
       % squares solution
% ss: step size for trimming
% tf: trimmed fraction (estimate of how many data points are corrupted)

[m,n] = size(A);
lambda = kappa*tau^2/2;

h = round((1-tf)*m);
iter = 0;
err = 10;

w_old = A\dxdt;
v_old = ones(m,1);
% -- Solve SINDy regression problem using SR3
while err > 10^(-8) && iter < maxiter
    iter = iter + 1;
    Theta_c = A.*repmat(v_old,1,n);
    H = Theta_c.'*A+ kappa *eye(n); 
    % -- update coefficients
    coeffs = H\(Theta_c.'*dxdt + kappa*eye(n)*w_old);
    % -- update w by taking l_0 prox 
    w_new = coeffs;
    w_new(abs(w_new)<=tau) = 0;
    % -- update v with gradient step and project onto capped h-simplex
    R2 = 0.5*sum((dxdt - A*coeffs).^2,2);
    v_new = proj_csimplex(v_old - ss*R2, h);
    
    % -- compute error
    err = norm(w_new - w_old)*kappa +norm(v_new-v_old)/ss;
    w_old = w_new; 
    v_old = v_new;
end

Xi_full = coeffs;
Xi_sparse = w_old;
trimmed_points = v_old;
end