% Exercise 3, analysis of Nelder-Mead method 
% Script for Evaluating the Nelder-Mead method.
% This script implements and tests the Nelder-Mead algorithm
% on three multidimensional benchmark functions:
%  - Chained Rosenbrock's Function
%  - Wood's Function
%  - Powell's Function
%
% For each function, simulations are performed to analyze the convergence, 
% computational cost, and empirical rates of convergence of the method in
% various dimensions (10,25,50) and multiple initial points (one suggested from 
% PDF and 10 points from the ipercube).
% It exploit NelderMead_for10Points() function

%% Chained Rosenbrook per dim = 10 %%
clc
clear all
close all
rng(min(343341,343428))  % ho diviso in sezioni quindi mettere a posto la questione seed 
dim = 3;
disp("dimensione:")
disp(dim)
f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

x1_opt = ones(1,dim);

[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f1_ros,x1_rosenbrock,x1_opt);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)  %fare comunque una media perchè vengono diversi a ogni iterazione nonostante il seed
disp("Empirical rate of convergence")
% for i = 1:length(lista_rates2)
%     disp(['Column ', num2str(i), ':']);
%     disp(lista_rates2{i});  % Mostra il contenuto della cella i
% end
% for i = 1:length(lista_rates1)
%     disp(['Column ', num2str(i), ':']);
%     disp(lista_rates1{i});  % Mostra il contenuto della cella i
% end
% vedi inizio di pag 26 di della santa
%% Chained Rosenbrook per dim = 25 %%
clc
clear all
close all
dim = 25;
disp("dimensione:")
disp(dim)
f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

x1_opt = ones(1,dim);

[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f1_ros,x1_rosenbrock,x1_opt);
% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)
for i = 1:length(lista_rates2)
    disp(['Column ', num2str(i), ':']);
    disp(lista_rates2{i});  % Mostra il contenuto della cella i
end
%% Chained Rosenbrook per dim = 50 %%
clc
clear all
close all
dim = 50;
disp("dimensione:")
disp(dim)
f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

x1_opt = ones(1,dim);

[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f1_ros,x1_rosenbrock,x1_opt);
% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)
for i = 1:length(lista_rates2)
    disp(['Column ', num2str(i), ':']);
    disp(lista_rates2{i});  % Mostra il contenuto della cella i
end
%% Wood function per dim = 10 %%
clc
clear all
close all
dim = 5;
disp("dimensione:")
disp(dim)

f2_wood = @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x2_wood = zeros(1,dim); 
x2_wood(mod(n,2) == 1 & n <= 4) = -3;
x2_wood(mod(n,2) == 1 & n > 4) = -2; 
x2_wood(mod(n,2) == 0 & n <= 4) = -1;
x2_wood(mod(n,2) == 0 & n > 4) = 0;   
x2_opt = ones(1,dim);
disp(x2_wood)
[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f2_wood,x2_wood,x2_opt);

% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

%% Wood function per dim = 25 %%
clc
clear all
close all
dim = 25;
disp("dimensione:")
disp(dim)

f2_wood = @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x2_wood = zeros(1,dim); 
x2_wood(mod(n,2) == 1 & n <= 4) = -3;
x2_wood(mod(n,2) == 1 & n > 4) = -2; 
x2_wood(mod(n,2) == 0 & n <= 4) = -1;
x2_wood(mod(n,2) == 0 & n > 4) = 0;   
x2_opt = ones(1,dim);
disp(x2_wood)

[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f2_wood,x2_wood,x2_opt);

% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)
for i = 1:length(lista_rates2)
    disp(['Column ', num2str(i), ':']);
    disp(lista_rates2{i});  % Mostra il contenuto della cella i
end
%% Wood function per dim = 50 %%
clc
clear all
close all
dim = 50;
disp("dimensione:")
disp(dim) 

f2_wood = @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x2_wood = zeros(1,dim); 
x2_wood(mod(n,2) == 1 & n <= 4) = -3;
x2_wood(mod(n,2) == 1 & n > 4) = -2; 
x2_wood(mod(n,2) == 0 & n <= 4) = -1;
x2_wood(mod(n,2) == 0 & n > 4) = 0;   
x2_opt = ones(1,dim);
disp(x2_wood)

[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f2_wood,x2_wood,x2_opt);

% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

%% Powell function per dim = 10 %%
clc
clear all
close all
dim = 4; % non minore di 4
disp("dimensione:")
disp(dim)

f3_powell = @(x) sum(arrayfun(@(j) ...
    (x(2*j-1) + 10*x(2*j))^2 + 5*(x(2*j+1) - x(2*j+2))^2 + ...
    (x(2*j) - 2*x(2*j+1))^4 + 10*(x(2*j-1) - x(2*j+2))^4, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x3_powell = zeros(1, dim); 
x3_powell(mod(n,4) == 1) = 3;   
x3_powell(mod(n,4) == 2) = -1; 
x3_powell(mod(n,4) == 3) = 0;  
x3_powell(mod(n,4) == 0) = 1;  
x3_opt = zeros(1,dim);

[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f3_powell,x3_powell,x3_opt);

% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

%% Powell function per dim = 25 %%
clc
clear all
close all
dim = 25;
disp("dimensione:")
disp(dim)

f3_powell = @(x) sum(arrayfun(@(j) ...
    (x(2*j-1) + 10*x(2*j))^2 + 5*(x(2*j+1) - x(2*j+2))^2 + ...
    (x(2*j) - 2*x(2*j+1))^4 + 10*(x(2*j-1) - x(2*j+2))^4, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x3_powell = zeros(1, dim); 
x3_powell(mod(n,4) == 1) = 3;   
x3_powell(mod(n,4) == 2) = -1; 
x3_powell(mod(n,4) == 3) = 0;  
x3_powell(mod(n,4) == 0) = 1;  
x3_opt = zeros(1,dim);

[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f3_powell,x3_powell,x3_opt);

% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

%% Powell function per dim = 50 %%
clc
clear all
close all
dim = 50;
disp("dimensione:")
disp(dim)

f3_powell = @(x) sum(arrayfun(@(j) ...
    (x(2*j-1) + 10*x(2*j))^2 + 5*(x(2*j+1) - x(2*j+2))^2 + ...
    (x(2*j) - 2*x(2*j+1))^4 + 10*(x(2*j-1) - x(2*j+2))^4, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x3_powell = zeros(1, dim); 
x3_powell(mod(n,4) == 1) = 3;   
x3_powell(mod(n,4) == 2) = -1; 
x3_powell(mod(n,4) == 3) = 0;  
x3_powell(mod(n,4) == 0) = 1;  
x3_opt = zeros(1,dim);

[vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f3_powell,x3_powell,x3_opt);

% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)