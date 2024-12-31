% Function: ComparePar
% This function evaluates the performance of the Nelder-Mead optimization 
% method using different parameter values (rho, sigma, gamma, chi) to 
% identify the best configuration for two goals:
% 1. Minimizing the number of iterations required for convergence.
% 2. Minimizing the error relative to the optimal solution.
% In each configuration a parameter is fixed.
%
% Inputs:
% dim           - Dimension of the optimization problem.
% f             - Objective function to be minimized.
% initial_point - Initial point for the Nelder-Mead method.
% x_opt         - Known optimal solution (used to compute distance).
% l             - Number of parameter configurations to test.
% rho           - Fixed reflection parameter.
% rho_vec       - Vector of reflection parameter values to test.
% sigma         - Fixed contraction parameter.
% sigma_vec     - Vector of contraction parameter values to test.
% gamma         - Fixed expansion parameter.
% gamma_vec     - Vector of expansion parameter values to test.
% chi           - Fixed external contraction parameter.
% chi_vec       - Vector of external contraction parameter values to test.
%
% Outputs:
% best_par_for_conv: vector with the value of each parameter which lead to
% the better convergence point compare with the optimal point
% best_par_for_conv = [rho_best, sigma_best, gamma_best, chi_best] 
%

function [best_par_for_conv,best_par_for_time] = comparePar(dim, f, initial_point,x_opt, l, rho,rho_vec,sigma,sigma_vec,gamma, gamma_vec,chi,chi_vec)
best_par_for_conv = [0,0,0,0];
best_par_for_time = [0,0,0,0];

num_par = 4;
kmax = 10000;
tol_simplex = 1e-07; 
tol_varf = 1e-07;  
matrix_value_par = [rho_vec; 
                    sigma_vec;
                    gamma_vec;
                    chi_vec];

[simplex_initial, flag2] = NelderMead_simplex(dim, initial_point);

k_vec_xPar = zeros(num_par,l);
dist_from_opt = zeros(num_par,l);

for i=1:l
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho_vec(i), chi, gamma, dim, sigma, tol_simplex, tol_varf);
    k_vec_xPar(1,i) = k;
    dist_from_opt(1,i) = norm(x_bar(end,:) - x_opt);
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi, gamma, dim, sigma_vec(i), tol_simplex, tol_varf);
    k_vec_xPar(2,i) = k;
    dist_from_opt(2,i) = norm(x_bar(end,:) - x_opt);
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi, gamma_vec(i), dim, sigma, tol_simplex, tol_varf);
    k_vec_xPar(3,i) = k;
    dist_from_opt(3,i) = norm(x_bar(end,:) - x_opt);
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi_vec(i), gamma, dim, sigma, tol_simplex, tol_varf);
    k_vec_xPar(4,i) = k;
    dist_from_opt(4,i) = norm(x_bar(end,:) - x_opt);
end
disp("k")
disp(k_vec_xPar)
disp("conv")
disp(dist_from_opt)

% find value of par that minimise the number of iteration
for i = 1:num_par
    disp(i)
    [minimo, idx_conv] = min(k_vec_xPar(i,:)); 
    disp([minimo, idx_conv])
    best_par_for_time(i) = matrix_value_par(i, idx_conv);  
end

% find value of par that minimise the error of convergence
for i = 1:num_par
    [minimo, idx_time] = min(dist_from_opt(i,:));  
    best_par_for_conv(i) = matrix_value_par(i, idx_time);  
end

