function [k1, x_bar1, k_10_points,x_bar_10_points] = NelderMead_for_10Points(dim,f,x_initial)
% vedere se mettere i parametri per nelder mead in imput

% creazione 10 punti dall'ipercubo
l_bound = x_initial - 1;
u_bound = x_initial + 1; 
M_ten_initial_points = zeros(dim, 10); % contains ten points in each col
for i=1:dim
    coord_random = l_bound(i)*ones(1,10) + (u_bound(i) - l_bound(i)) * rand(1, 10);
    M_ten_initial_points(i,:) = coord_random;
end

% Parametri per la funzione
kmax = 100000;
rho = 1;
sigma = 1/2;
gamma = 1/2;
chi = 2;
tol_simplex = 1e-06;  % tolleranza su simplesso vedere che 10 alla meno 7
tol_varf = 1e-06;   % tolleranza su f

% Nelder Mead method with x initial points
[simplex_initial, flag] = NelderMead_simplex(dim, x_initial);
[k1, simplex,x_bar1, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi, gamma, dim, sigma, tol_simplex, tol_varf);


% Nelder Mead method with 10 initial points from ipercube
x_bar_10_points = zeros(10,dim + 1);  %inserisco il flag a fine riga
k_10_points = zeros(1,10);
for i = 1:10
    initial_point = M_ten_initial_points(:,i)';
    disp(" ")
    disp(i)
    [simplex_initial, flag] = NelderMead_simplex(dim, initial_point);
    [k, simplex,x_bar, flag]  = nelder_mead(f, simplex_initial, kmax, rho, chi, gamma, dim, sigma, tol_simplex, tol_varf);
    row_x_bar_10_points = [x_bar(end,:), flag];
    x_bar_10_points(i,:) = row_x_bar_10_points ;
    k_10_points(1,i) = k;
    % distance_from_opt = norm(initial_point-ones(1,dim)); non sembra centrare
end
end
