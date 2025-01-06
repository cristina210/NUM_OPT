function [simplex_initial, flag] = NelderMead_simplex(dim, initial_point)
flag = 0;
simplex_initial = zeros(dim+1, dim); 
for i = 1:(dim+1)
    for j = 1:(dim)
        if i >= 2 && i == j + 1
            simplex_initial(i,j) = 1;
        else
            simplex_initial(1,:) = initial_point;
        end
    end
end
% Controllo che sia un simplesso non degenere (differenza tra un un punto e il primo
% punto sono lin. ind.)
vettore_differenze = zeros(dim,dim);
for i = 2:dim+1
    vettore_differenze(i-1,:) = simplex_initial(1,:) - simplex_initial(i,:);
end
if rank(vettore_differenze) ~= dim
    disp("Il simplesso iniziale è degenere")
    flag = 1;
end
end

