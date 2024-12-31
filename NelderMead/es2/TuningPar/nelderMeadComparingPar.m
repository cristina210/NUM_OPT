clc
clear all
close all

% Each parameter is tested with 7 different values, while others remain
% constant. Then number of iteration before convergence and distance of the
% final barycenter from the optimal solution is compared.

dim = 2;  
initial_point = [-1,1.2];  
% initial_point = [1.2,1.2];
x_opt = [1,1];
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ; % Rosenbrock function

num_par = 4;
l = 5;  % number of possible value of each parameter to be tested
kmax = 10000;
rho = 1; 
rho_vec = [0.2,0.5,1,1.5,1.7];
sigma = 1/2;
sigma_vec = [0.1,0.3,0.5,0.7,0.9];
gamma = 1/2;
gamma_vec = [0.1,0.3,0.5,0.7,0.9];
chi = 2;
chi_vec = [1, 1.5, 2, 2.5, 3];

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
disp(k_vec_xPar)
disp(dist_from_opt)

for j=1:num_par
 figure;
 plot(matrix_value_par(j,:),k_vec_xPar(j,:),'-go','MarkerSize',9,'MarkerFaceColor','g')
 hold on
 plot(matrix_value_par(j,4),k_vec_xPar(j,4), '--*','Color','k')
 title('Number of Iterations');  
 if j == 1
        xlabel('rho Value');
    elseif j == 2
        xlabel('sigma Value');
    elseif j == 3
        xlabel('gamma Value');
    elseif j == 4
        xlabel('chi Value');
    end
 ylabel('Number of Iterations');
 figure;
 plot(matrix_value_par(j,:),dist_from_opt(j,:),'-ro','MarkerSize',9,'MarkerFaceColor','r')
 hold on
 plot(matrix_value_par(j,4),dist_from_opt(j,4), '--*','Color','k')
 title('Distance from Optimal Solution')
 if j == 1
        xlabel('rho Value');
    elseif j == 2
        xlabel('sigma Value');
    elseif j == 3
        xlabel('gamma Value');
    elseif j == 4
        xlabel('chi Value');
    end
  ylabel('Number of Iterations');
end