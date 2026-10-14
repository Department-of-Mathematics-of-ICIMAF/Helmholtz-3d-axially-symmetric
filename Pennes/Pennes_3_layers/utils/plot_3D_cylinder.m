function plot_3D_cylinder(Rho, Z, U, theta)
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
surf(X_cyl, Y_cyl, Z_cyl, U_cyl,'EdgeColor', 'interp'); %'FaceColor','interp' , 
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Superficie rotada');

%% Graficar las caras del volumen
plot_face_3D(Rho, Z, U, theta, true);

%% Graficar el Eje Rotacion
[X_line, Y_line] = meshgrid([0 0], [0 0]);
min_z = min(min(Z));
max_z = max(max(Z));
Z_line = [min_z min_z;max_z*1.2 max_z*1.2];
hold on
view([48,29])
surf(X_line, Y_line, Z_line, zeros(2,2), 'LineWidth', 1);
% min_value = min(U);
% if length(min_value)>1
%     min_value = min(min_value);
% end
% max_value = max(U);
% if length(max_value)>1
%     max_value = max(max_value);
% end
% colorbar('Ticks',linspace(min_value,max_value,6));

hold off

end