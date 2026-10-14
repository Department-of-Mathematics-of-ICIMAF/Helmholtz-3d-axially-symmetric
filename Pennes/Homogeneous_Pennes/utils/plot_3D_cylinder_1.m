function plot_3D_cylinder_1(Rho, Z, U, theta, min_value, max_value, title_bar, z1, z2, z3)
% Construir grafico por rotacion del plano Rho-Z desde el angulo 0 hasta el
% angulo theta
% ENTRADA:
% Rho - Matriz de (npoints, npoints) que contiene los radios respecto al
% eje Z
% Z - Matriz de (npoints, npoints) que contiene los valores respecto al eje Z
% U - Matriz de (npoints, npoints) que contiene los valores de la funcion
% en el punto (rho, z)
% theta - angulo final de la rotacion



%% Extraer perfiles de Rho, de Z y de U
[n,m] = size(Rho);
rho_perfil = Rho(:,1);
z_perfil = Z(:,1);
u_perfil = U(:,1);
for i=2:n-1
    rho_perfil = [rho_perfil; Rho(n,i)];
    z_perfil = [z_perfil; Z(n, i)];
    u_perfil = [u_perfil; U(n, i)];
end
rho_perfil = [rho_perfil; Rho(n:-1:1,n)];
z_perfil = [z_perfil; Z(n:-1:1,n)];
u_perfil = [u_perfil; U(n:-1:1,n)];

%% Usar cylinder para crear la superficie
[X_cyl, Y_cyl, Z_cyl, U_cyl] = cylinder_xyz(rho_perfil, 100, z_perfil, u_perfil, theta);

% figure;
surf(X_cyl, Y_cyl, Z_cyl, U_cyl,'EdgeColor', 'interp', 'FaceAlpha', 1); %'FaceColor','interp' , 
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

%% Graficar las caras del volumen
plot_face_3D_1(Rho, Z, U, theta, true, z1, z2, z3);
% colorbar('Ticks',linspace(min_value,max_value,6));

%% Graficar el Eje Rotacion
hold on
min_z = min(min(Z));
max_z = max(max(Z));
plot3([0 0], [0 0], [min_z max_z*1.2], 'LineStyle', '-', 'LineWidth', 2, 'Color', 'k'); % 
hold off

%% Graficar los limites de capas
% hold on
% max_x = max(max(Rho));
% min_x = min(min(Rho));
% plot3([min_x max_x], [0 0], [z1 z1], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'm'); % 
% plot3([min_x max_x], [0 0], [z2 z2], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'm'); % 
% % plot3([0 max_rho], [0 0], [z3 z3], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'm'); % 
% hold off

cb = colorbar;
cb.Limits = [min_value, max_value];
% cb.Label.String = title_bar;
view([48,29])
end