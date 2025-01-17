function [vec_e_k] = stagnation_func(x_bar)
vec_e_k = zeros(1,length(x_bar)-1);
for i=2:length(x_bar)
    e_k = norm(x_bar(i,:) - x_bar(i-1,:));
    vec_e_k(1,i-1) = e_k;
end
figure;
plot(1:length(vec_e_k), vec_e_k, 'g', 'LineWidth', 2)
title('Study of  stagnation point','FontSize', 16)
xlabel('iteration', 'FontSize', 14)
ylabel('||x(k+1) - x(k)||', 'FontSize', 14)

