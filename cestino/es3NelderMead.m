clc
clear all
close all

%% Rosenbrook Chained per dim = 10 %%
clc
clear all
close all
rng(min(343341,343428)) 
dim = 10;
disp("dimensione:")
disp(dim)
f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

[k1, x_bar1, k_10_points,x_bar_10_points] = NelderMead_for_10Points(dim,f1_ros,x1_rosenbrock);
% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1(end,:))
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
%% Rosenbrook Chained per dim = 25 %%
clc
clear all
close all
dim = 25;
disp("dimensione:")
disp(dim)
f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

[k1, x_bar1, k_10_points,x_bar_10_points] = NelderMead_for_10Points(dim,f1_ros,x1_rosenbrock);
% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1(end,:))
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
%% Rosenbrook Chained per dim = 50 %%
clc
clear all
close all
dim = 50;
disp("dimensione:")
disp(dim)
f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

[k1, x_bar1, k_10_points,x_bar_10_points] = NelderMead_for_10Points(dim,f1_ros,x1_rosenbrock);
% i parametri per il nelder mead li ho messi dentro la funzione
disp("Convergence point from initial point suggested from pdf file:")
disp(x_bar1(end,:))
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_bar_10_points)
%% Wood function per dim = 10 %%
clc
clear all
close all
dim = 10;
disp("dimensione:")
disp(dim)

f2_wood = @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

x2_wood = arrayfun(@(i) ...
    (-3*(mod(i,2)==1 && i>4) - 1*(mod(i,2)==1 && i<=4)) + ...
    (0*(mod(i,2)==0 && i>4) - 1*(mod(i,2)==0 && i<=4)), 1:n);

% to be continue
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

x2_wood = arrayfun(@(i) ...
    (-3*(mod(i,2)==1 && i>4) - 1*(mod(i,2)==1 && i<=4)) + ...
    (0*(mod(i,2)==0 && i>4) - 1*(mod(i,2)==0 && i<=4)), 1:n);

% to be continue
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

x2_wood = arrayfun(@(i) ...
    (-3*(mod(i,2)==1 && i>4) - 1*(mod(i,2)==1 && i<=4)) + ...
    (0*(mod(i,2)==0 && i>4) - 1*(mod(i,2)==0 && i<=4)), 1:n);

% to be continue
%% Powell function per dim = 10 %%
clc
clear all
close all
dim = 10;
disp("dimensione:")
disp(dim)

f3_powell = @(x) sum(arrayfun(@(j) ...
    (x(2*j-1) + 10*x(2*j))^2 + 5*(x(2*j+1) - x(2*j+2))^2 + ...
    (x(2*j) - 2*x(2*j+1))^4 + 10*(x(2*j-1) - x(2*j+2))^4, ...
    1:(length(x)-2)/2));

x3_powell = arrayfun(@(i) ...
    (3*(mod(i,4)==1) - 1*(mod(i,4)==2) + ...
     0*(mod(i,4)==3) + 1*(mod(i,4)==0)), 1:n);

% to be continue
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

x3_powell = arrayfun(@(i) ...
    (3*(mod(i,4)==1) - 1*(mod(i,4)==2) + ...
     0*(mod(i,4)==3) + 1*(mod(i,4)==0)), 1:n);

% to be continue
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

x3_powell = arrayfun(@(i) ...
    (3*(mod(i,4)==1) - 1*(mod(i,4)==2) + ...
     0*(mod(i,4)==3) + 1*(mod(i,4)==0)), 1:n);

% to be continue