function nelderMead_picture2D(f, x_interval, y_interval, x_bar1, x_bar2, initial_point1, initial_point2)
[X, Y] = meshgrid(x_interval, y_interval);   
Z = f(X, Y);
figure;
imagesc(x_interval, y_interval, Z);
set(gca, 'YDir', 'normal'); 
colorbar;
set(gca, 'ColorScale', 'log');  
hold on;
contour(X, Y, Z, 50, 'LineColor', 'k');  
plot(initial_point1(1), initial_point1(2), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); 
plot(initial_point2(1), initial_point2(2), 'cp', 'MarkerSize', 12, 'MarkerFaceColor', 'c'); 
plot(x_bar1(:, 1), x_bar1(:, 2), 'ro-', 'LineWidth', 1.3, 'MarkerSize', 3.5);
plot(x_bar2(:, 1), x_bar2(:, 2), 'co-', 'LineWidth', 1.3, 'MarkerSize', 3.5);
plot(1, 1, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 5); 
hold off;
xlabel('x');
ylabel('y');
title('Funzione di Rosenbrock - Mappa di Calore e Contorni');
end

