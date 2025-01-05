% Function: NelderMead_for_10Points
% Runs the Nelder-Mead method for a given function and  
% for a fixed initial point and 10 random points within a hypercube build 
% from the fixed initial point. 
% This function measures performance, including convergence rates and execution times, of
% the algorithm.
%
% Inputs:
% - dim: The dimension of the problem.
% - f: The objective function to minimize
% - x_initial: The initial point 
% - x_opt: The known optimal solution for the objective function
%
% Outputs:
% - vec_time: A 1x11 vector of elapsed times for the Nelder-Mead algorithm.
%     - vec_time(1): time for the single initial point.
%     - vec_time(2:end): times for each of the 10 random initial points.
% - k1: Number of iterations required for convergence using the single initial point.
% - x_bar1: Final solution (final barycenter of the simplex) for the single initial point.
% - k_10_points: A 1x10 vector containing the number of iterations required for convergence 
%                  for each of the 10 random initial points.
% - x_bar_10_points: A 10x(dim+1) matrix where each row contains:
%     - The final solution (final barycenter of the simplex) for each random initial point.
%     - A success/failure flag as the last element of the row.
% - lista_rates1: Cell array containing convergence rate vectors for each random point 
%                   using an alternative computation method (currently unused in the code).
% - lista_rates2: Cell array containing convergence rate vectors for each random point 
%                   based on the distance to the known optimal solution.
% This function makes use of some functions such as compute_exp_rate_conv2()

function [vec_time, k1, x_bar1, k_10_points,x_bar_10_points, lista_rates1, lista_rates2] = NelderMead_for_10Points(dim,f,x_initial,x_opt)

% Parameters for Nelder Mead
kmax = 100000;
rho = 1;
sigma = 1/2;
gamma = 1/2;
chi = 2;
tol_simplex = 1e-07;  
tol_varf = 1e-07; 
vec_time = zeros(1,11); % contains computational costs

% Nelder Mead method with suggested initial points
[simplex_initial, flag] = NelderMead_simplex(dim, x_initial);
tic
[k1, simplex,x_bar1, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi, gamma, dim, sigma, tol_simplex, tol_varf);
vec_time(1,1) = toc;
%figure;
% compute_exp_rate_conv(x_bar1, k1);
figure;
vec_rate2 = compute_exp_rate_conv2(x_bar1, k1, x_opt);
x_bar1 = x_bar1(end,:);

% Nelder Mead method with 10 random initial points from ipercube
% Generate points
l_bound = x_initial - 1;
u_bound = x_initial + 1; 
M_ten_initial_points = zeros(dim, 10); 
for i=1:dim
    coord_random = l_bound(i)*ones(1,10) + (u_bound(i) - l_bound(i)) * rand(1, 10);
    M_ten_initial_points(i,:) = coord_random;
end
lista_rates1 = {};
lista_rates2 = {};

% Nelder mead for each point
x_bar_10_points = zeros(10,dim + 1);  
k_10_points = zeros(1,10);

for i = 1:10
    i
    initial_point = M_ten_initial_points(:,i)';
    disp(" ")
    disp(i)
    [simplex_initial, flag] = NelderMead_simplex(dim, initial_point);
    tic
    [k, simplex, x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi, gamma, dim, sigma, tol_simplex, tol_varf);
    vec_time(1,i+1) = toc;
    row_x_bar_10_points = [x_bar(end,:), flag];
    x_bar_10_points(i,:) = row_x_bar_10_points ;
    k_10_points(1,i) = k;
    %figure;
    %vec_rate1_i = compute_exp_rate_conv(x_bar, k);
    %lista_rates1{end+1} = vec_rate1_i;
    figure;
    vec_rate2_i = compute_exp_rate_conv2(x_bar, k, x_opt);
    lista_rates2{end+1} = vec_rate2_i;
end
end

