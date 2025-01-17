% Exercise 2, Nelder Mead method
% In this script the Nelder Mead method is applied to minimize the Rosenbrock
% function in two dimensions with two different starting point. 

clc
clear all
close all

dim = 2;  
initial_point1 = [1.2,1.2];
initial_point2 = [-1,1.2];
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ; % Rosenbrock function
x_opt = [1,1];

% Tuned parameters for point 1
kmax = 10000;
rho1 = 0.5;
sigma1 = 0.1;
gamma1 = 0.3;
chi1 = 2.3;

% Tuned parameters for point 2
rho2 = 1;
sigma2 = 0.6;
gamma2 = 0.5;
chi2 = 2;

% tol for first stopping criterion
tol_simplex = 1e-07;  
% tol for second stopping criterion
tol_varf = 1e-07;   % 1e-17

% Create initial simplex
[simplex_initial1, flag1] = NelderMead_simplex(dim, initial_point1);
[simplex_initial2, flag2] = NelderMead_simplex(dim, initial_point2);

% Nelder mead
tic;
[k1, simplex1,x_bar1, flag1]  = nelder_mead(f, simplex_initial1, kmax, rho1, chi1, gamma1, dim, sigma1, tol_simplex, tol_varf);
time1 = toc;
tic;
[k2, simplex2,x_bar2, flag2]  = nelder_mead(f, simplex_initial2, kmax, rho2, chi2, gamma2, dim, sigma2, tol_simplex, tol_varf);
time2 = toc;
x_bar2(end,:);

% Output
disp("[1.2,1.2] initial point, barycenter of the convergence simplex")
disp(x_bar1(end,:))
disp("[-1,1.2] initial point, barycenter of the convergence simplex")
disp(x_bar2(end,:))
disp("[1.2,1.2] initial point, number iteration before convergence")
disp(k1)
disp("[-1,1.2] initial point, number iteration before convergence")
disp(k2)
disp("[1.2,1.2] initial point, computational costs")
disp(time1)
disp("[-1,1.2] initial point, computational costs")
disp(time2)

vec_rate1 = compute_exp_rate_conv2(x_bar1, k1, x_opt);   % order of convergence
vec_rate2 = compute_exp_rate_conv2(x_bar2, k2, x_opt);    % order of convergence
vec_increments1 = stagnation_func(x_bar1);
vec_increments2 = stagnation_func(x_bar2);

% Picture
f = @(x, y) 100*(y - x.^2).^2 + (1 - x).^2;
x_interval = linspace(-2, 2, 500);  
y_interval = linspace(-1, 3, 500); 
nelderMead_picture2D(f, x_interval, y_interval, x_bar1, x_bar2, initial_point1, initial_point2)

