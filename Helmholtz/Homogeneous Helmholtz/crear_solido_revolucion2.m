function crear_solido_revolucion2(ax, X, Y, Z, angulo_total, clim)
    % Agrega un sólido de revolución parcial (cuña) con sus caras planas laterales.
    % ax: ejes donde graficar
    % X, Y: matrices de coordenadas (radio vs altura) en [0,1]
    % Z: matriz de valores (mapa de colores)
    % angulo_total: ángulo de revolución en radianes (desde 0 hasta ese valor)
    % clim: (opcional) vector [min max] para escala de colores global

    if nargin < 6
        clim = [];
    end

    % Verificar tamaños
    if ~isequal(size(X), size(Y), size(Z))
        error('Las matrices X, Y y Z deben tener el mismo tamaño');
    end

    % Resolución angular de la superficie curva
    num_angulos = 100;
    angulos = linspace(0, angulo_total, num_angulos);

    % Datos del perfil
    radios = X(:);
    alturas = Y(:);
    valores = Z(:);

    % ---- Superficie curva (revolución) ----
    [Theta_grid, R_grid] = meshgrid(angulos, radios);
    [~, H_grid] = meshgrid(angulos, alturas);

    X_curva = R_grid .* cos(Theta_grid);
    Y_curva = R_grid .* sin(Theta_grid);
    Z_curva = H_grid;
    C_curva = repmat(valores, 1, num_angulos);

    surf(ax, X_curva, Y_curva, Z_curva, C_curva, 'EdgeColor', 'interp','FaceAlpha',1);
    hold(ax, 'on');

    % ---- Cara plana en ángulo 0 (θ = 0) ----
    % En θ=0: x = r, y = 0, z = h
    X0 = X;           % matriz igual que X
    Y0 = zeros(size(X));
    Z0 = Y;           % altura
    C0 = Z;           % color según el valor original

    surf(ax, X0, Y0, Z0, C0, 'EdgeColor', 'interp','FaceAlpha',1);

    % ---- Cara plana en ángulo final (θ = angulo_total) ----
    cos_ang = cos(angulo_total);
    sin_ang = sin(angulo_total);
    Xf = X * cos_ang;
    Yf = X * sin_ang;
    Zf = Y;           % misma altura
    Cf = Z;           % mismo color

    surf(ax, Xf, Yf, Zf, Cf, 'EdgeColor', 'interp','FaceAlpha',1);

    % ---- Escala de colores global, si se proporciona ----
    if ~isempty(clim)
        caxis(ax, clim);
    end
end