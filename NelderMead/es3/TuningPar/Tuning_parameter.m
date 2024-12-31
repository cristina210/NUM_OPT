% Tuning of parameters rho, sigma, gamma, chi used in Nelder Mead method.
% 3 different values of each parameters are compared in each dimensione =
% 10,25,50. In each section a different function of ex 3 is take into
% account together with the initial starting point suggested from the pdf
% file.
% The comparison between parameters is accomplished through the study of
% two quantities: number of iteration to converge and error of convergence
% (distance from optimal point)
% The choice for the parameter is done take into consideration these
% two quantities: best parameter shoud minimize both of them.
% It is exploited function "comparePar()"

%% Tuning in Chained Rosenbrock
clc
clear all
close all

% Value of parameters compared
l = 3; 
rho_vec = [0.25,1,1.75];
sigma_vec = [0.2,0.5,0.8];
gamma_vec = [0.2,0.5,0.8];
chi_vec = [1,2,3];

% Comparison in different dimension
vec_dim = [5];
% vec_dim = [10,25,50]

% Rosenbrock function
for j = 1:length(vec_dim)
    
    disp("Dimension:")
    dim = vec_dim(j);
    disp(dim)

    % Information of function and initial point (suggested from PDF)
    f= @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));
    initial_point = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);
    x_opt = ones(1,dim);

    % Nelder Mead
    [pos_bestRho, pos_bestSigma, pos_bestGamma, pos_bestChi] = comparePar2(dim, f, initial_point,x_opt, l, rho_vec,sigma_vec, gamma_vec,chi_vec);
    disp("Best configuration:")
    disp("value of rho")
    disp(rho_vec(pos_bestRho))
    disp("value of sigma")
    disp(sigma_vec(pos_bestSigma))
    disp("value of gamma")
    disp(gamma_vec(pos_bestGamma))
    disp("value of chi")
    disp(chi_vec(pos_bestChi))
    %Altra opzione se l'altra è troppo lenta)
    %[best_par_for_conv,best_par_for_time] = comparePar(dim, f, initial_point,x_opt, l, rho,rho_vec,sigma,sigma_vec,gamma, gamma_vec,chi,chi_vec);
    %disp(best_par_for_time)
    %disp(best_par_for_conv)
end  

%% Tuning in Wood's function
% stessa cosa per wood e powell
clc
clear all
close all

% Value of parameters compared
l = 3; 
rho = 1; 
rho_vec = [0.25,1,1.75];
sigma = 1/2;
sigma_vec = [0.2,0.5,0.8];
gamma = 1/2;
gamma_vec = [0.2,0.5,0.8];
chi = 2;
chi_vec = [1,2,3];

% Comparison in different dimension
vec_dim = [10,15];
% vec_dim = [10,25,50]

% Wood function
for j = 1:length(vec_dim)
    dim = vec_dim(j);

    % Information of function and initial point (suggested from PDF)
    f= @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

    n = 1:dim; 
    initial_point = zeros(1,dim); 
    initial_point(mod(n,2) == 1 & n <= 4) = -3;
    initial_point(mod(n,2) == 1 & n > 4) = -2; 
    initial_point(mod(n,2) == 0 & n <= 4) = -1;
    initial_point(mod(n,2) == 0 & n > 4) = 0;   
    x_opt = ones(1,dim);

    % stessa cosa
end 
