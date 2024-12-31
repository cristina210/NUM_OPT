%%%%%%%%%%%%%%%%%%%%%%%%%%% Nelder mead %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

% Dati da fornire
dim = 2;   % dimensione del dominio della funzione
initial_point1 = [1.2,1.2];
initial_point2 = [-1,1.2];
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ; % funzione di Rosenbrock

% Parametri per la funzione
kmax = 1000;
rho = 1;
sigma = 1/2;
gamma = 1/2;
chi = 2;
tol_simplex = 1e-07;  % tolleranza su simplesso
tol_varf = 1e-07;   % tolleranza su f

% Creazione dei due simplessi iniziali 
% formato da punti: punto iniziale e i punti che formano una base
% ortonormale in R^{dim + 1}

% Create initial simplex
[simplex_initial1, flag1] = NelderMead_simplex(dim, initial_point1);
[simplex_initial2, flag2] = NelderMead_simplex(dim, initial_point2);

% Applico la funzione nelder_mead
% k ultima iterazione, simplex1 simplesso finale, x_bar = vettore che
% contiene i baricentri nelle varie iterazioni
tic;
[k1, simplex1,x_bar1, flag1]  = nelder_mead(f, simplex_initial1, kmax, rho, chi, gamma, dim, sigma, tol_simplex, tol_varf);
size(x_bar1)
time1 = toc
k1;
flag1;
x_bar1(end,:);
tic;
[k2, simplex2,x_bar2, flag2]  = nelder_mead(f, simplex_initial2, kmax, rho, chi, gamma, dim, sigma, tol_simplex, tol_varf);
time2 = toc
k2;
flag2;
x_bar2(end,:);

% Picture
f = @(x, y) 100*(y - x.^2).^2 + (1 - x).^2;
x_interval = linspace(-2, 2, 500);  
y_interval = linspace(-1, 3, 500); 
nelderMead_picture2D(f, x_interval, y_interval, x_bar1, x_bar2, initial_point1, initial_point2)

% Comparing
vec_rate1 = compute_exp_rate_conv(x_bar1, k1);
vec_rate2 = compute_exp_rate_conv(x_bar2, k2);
vec_rate2 
