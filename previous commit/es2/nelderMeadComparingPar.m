clc
clear all
close all

% Dati da fornire
dim = 2;   % dimensione del dominio della funzione
initial_point = [-1,1.2];   %viene meglio con questo
% initial_point = [1.2,1.2];
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ; % funzione di Rosenbrock

% Parametri per la funzione
l = 7;
num_par = 4;
kmax = 300;
rho = 1;
rho_vec = [0.4,0.5,0.7,1,1.3,1.5,1.7];
sigma = 1/2;
sigma_vec = [0.1,0.2,0.3,0.5,0.65,0.8,0.95];
gamma = 1/2;
gamma_vec = [0.2,0.3,0.4,0.5,0.65,0.8,0.95];
chi = 2;
chi_vec = [1,1.3, 1.5,2,2.3 2.5,3];
tol_simplex = 1e-07;  % tolleranza su simplesso
tol_varf = 1e-07;   % tolleranza su f
matrix_value_par = [rho_vec; 
                    sigma_vec;
                    gamma_vec;
                    chi_vec];

[simplex_initial, flag2] = NelderMead_simplex(dim, initial_point);

k_vec_xPar = zeros(num_par,l);
% row1: # iteration of each value of rho
% row2: # iteration of each value of sigma
% row3: # iteration of each value of gamma
% row4: # iteration of each value of rho

for i=1:l
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho_vec(i), chi, gamma, dim, sigma, tol_simplex, tol_varf);
    k_vec_xPar(1,i) = k;
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi, gamma, dim, sigma_vec(i), tol_simplex, tol_varf);
    k_vec_xPar(2,i) = k;
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi, gamma_vec(i), dim, sigma, tol_simplex, tol_varf);
    k_vec_xPar(3,i) = k;
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi_vec(i), gamma, dim, sigma, tol_simplex, tol_varf);
    k_vec_xPar(4,i) = k;
end
disp(k_vec_xPar)

for j=1:num_par
 figure;
 plot(matrix_value_par(j,:),k_vec_xPar(j,:),'-bo','MarkerSize',6,'MarkerFaceColor','b')
 hold on
 plot(matrix_value_par(j,4),k_vec_xPar(j,4), '--*','Color','r')
end