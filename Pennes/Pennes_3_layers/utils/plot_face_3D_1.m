function plot_face_3D_1(Rho, Z, U, theta_end, add_plot, z1, z2, z3)
% Crear m?ltiples cortes rotados

theta = [0 theta_end];


if add_plot ~= true
    figure;
    axis equal;
    colorbar;
    view([48,29])
    % title('M?ltiples cortes radiales');
else
    
end
hold on;
for k = 1:length(theta)
    % Coordenadas del slice rotado
    X_slice = Rho * cos(theta(k));
    Y_slice = Rho * sin(theta(k));
    Z_slice = Z;
    
%     surf(X_slice, Y_slice, Z_slice, U, 'EdgeColor', 'interp');
surf(X_slice, Y_slice, Z_slice, U, 'EdgeColor', 'interp', 'FaceAlpha', 0.5);
 
end

%% Graficar los limites de capas
max_x = max(max(X_slice));
min_x = min(min(X_slice));
max_y = max(max(Y_slice));
min_y = min(min(Y_slice));
hold on
plot3([min_x max_x], [min_y max_y], [z1 z1], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'm'); % 
plot3([min_x max_x], [min_y max_y], [z2 z2], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'm'); % 
% plot3([0 max_x], [0 max_y], [z3 z3], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'm'); % 
hold off
hold off
end