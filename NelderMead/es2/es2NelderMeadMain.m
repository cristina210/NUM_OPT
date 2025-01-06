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

% Parameters 
kmax = 1000;
rho = 1;
sigma = 1/2;
gamma = 1/2;
chi = 2;
tol_simplex = 1e-07;  
tol_varf = 1e-07;

% Create initial simplex
[simplex_initial1, flag1] = NelderMead_simplex(dim, initial_point1);
[simplex_initial2, flag2] = NelderMead_simplex(dim, initial_point2);

% Nelder mead
tic;
[k1, simplex1,x_bar1, flag1]  = nelder_mead(f, simplex_initial1, kmax, rho, chi, gamma, dim, sigma, tol_simplex, tol_varf);
size(x_bar1)
time1 = toc;
tic;
[k2, simplex2,x_bar2, flag2]  = nelder_mead(f, simplex_initial2, kmax, rho, chi, gamma, dim, sigma, tol_simplex, tol_varf);
time2 = toc;
x_bar2(end,:);

% Output
disp("Convergence point from first initial point:")
disp(x_bar1(end,:))
disp("Convergence point from second initial point:")
disp(x_bar2(end,:))
disp("Number iteration before convergence from first initial point:")
disp(k1)
disp("Number iteration before convergence from second initial point:")
disp(k2)
vec_rate1 = compute_exp_rate_conv2(x_bar1, k1, x_opt);   % order of convergence
vec_rate2 = compute_exp_rate_conv2(x_bar2, k2, x_opt);    % order of convergence

% Picture
f = @(x, y) 100*(y - x.^2).^2 + (1 - x).^2;
x_interval = linspace(-2, 2, 500);  
y_interval = linspace(-1, 3, 500); 
nelderMead_picture2D(f, x_interval, y_interval, x_bar1, x_bar2, initial_point1, initial_point2)

% COSA MANCA: tempi di convergenza studiarli in un altro file (prendere i
% tempi medi!!!)