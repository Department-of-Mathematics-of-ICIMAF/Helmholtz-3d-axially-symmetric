function plot_face_3D(Rho, Z, U, theta_end, add_plot)
% Crear m?ltiples cortes rotados

theta = [0 theta_end];


if add_plot ~= true
    figure;
    axis equal;
%     colorbar;
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
    
    surf(X_slice, Y_slice, Z_slice, U, 'EdgeColor', 'interp');
end
hold off
end