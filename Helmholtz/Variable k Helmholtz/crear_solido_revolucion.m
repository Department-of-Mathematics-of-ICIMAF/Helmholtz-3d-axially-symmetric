function crear_solido_revolucion(X, Y, Z)
    % Función para crear un sólido de revolución rotando el mapa de colores alrededor de x=0
    % X, Y: matrices de coordenadas (deben estar en el rango [0,1])
    % Z: matriz de valores (mapa de colores)
    
    % Verificar que las matrices tengan el mismo tamaño
    if ~isequal(size(X), size(Y), size(Z))
        error('Las matrices X, Y y Z deben tener el mismo tamaño');
    end
    
    % Número de puntos de rotación (resolución angular)
    num_angulos = 100;
    
    % Crear ángulos de rotación (rotación completa: 0 a 2π)
    angulos = linspace(0, 2*pi, num_angulos);
    
    % Preparar matrices para el sólido de revolución
    radios = X(:);  % Vector columna de radios
    alturas = Y(:); % Vector columna de alturas
    valores = Z(:); % Vector columna de valores Z
    
    % Crear malla de coordenadas cilíndricas
    [Theta_grid, R_grid] = meshgrid(angulos, radios);
    [~, H_grid] = meshgrid(angulos, alturas);
    
    % Convertir a coordenadas cartesianas 3D
    X_3D = R_grid .* cos(Theta_grid);
    Y_3D = R_grid .* sin(Theta_grid);
    Z_3D = H_grid;
    
    % Crear matriz de valores para la superficie
    Valores_grid = repmat(valores, 1, num_angulos);
    
    % Crear figura
    figure('Position', [100, 100, 1200, 500]);
    
    % Subplot 1: Mapa de colores original
    subplot(1, 2, 1);
    surf(X, Y, zeros(size(Z)), Z, 'EdgeColor', 'none');
    view(2);
    axis equal;
    xlabel('X (Radio)');
    ylabel('Y (Altura)');
    title('Mapa de colores original');
    colorbar;
    grid on;
    
    % Subplot 2: Sólido de revolución
    subplot(1, 2, 2);
    surf(X_3D, Y_3D, Z_3D, Valores_grid, 'EdgeColor', 'none');
    
    % Configurar vista 3D
    axis equal;
    axis vis3d;
    view(3);
    xlabel('X');
    ylabel('Y');
    zlabel('Z (Altura)');
    title('Sólido de revolución alrededor de x=0');
    colorbar;
    grid on;
    
    % Añadir iluminación
    camlight left;
    lighting gouraud;
    material dull;
    
    % Mostrar información
    fprintf('Sólido de revolución creado:\n');
    fprintf('- Puntos en el perfil: %d\n', length(radios));
    fprintf('- Resolución angular: %d puntos\n', num_angulos);
    fprintf('- Total de puntos: %d\n', numel(X_3D));
end