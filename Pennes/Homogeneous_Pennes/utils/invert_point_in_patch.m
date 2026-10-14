function [xi, eta, patch, info] = invert_point_in_patch(geometry, x, y, xi0, eta0)
    % geometry: la estructura geometry correspondiente
    % x, y: coordenadas del punto f?sico objetivo
    % xi0, eta0: estimaci?n inicial (por defecto (0.5, 0.5))
    if nargin < 4
        xi0 = 0.5;
        eta0 = 0.5;
    end
    
    npatch = numel(geometry);
    
    all_patches = [];
    all_xi      = [];
    all_eta     = [];
    all_res     = [];
    
    for j = 1:npatch
        geom_j = geometry(j);
        
        % Comprobar que el parche tiene los manejadores necesarios
        if ~isfield(geom_j, 'map') || ~isfield(geom_j, 'map_der')
            warning('El parche %d no tiene .map o .map_der, se omite.', j);
            continue;
        end
        
        [xi_j, eta_j, conv, res_j] = invert_point(geom_j, x, y, xi0, eta0);
        
        if conv
            all_patches(end+1) = j;      %#ok<AGROW>
            all_xi(end+1)      = xi_j;   %#ok<AGROW>
            all_eta(end+1)     = eta_j;  %#ok<AGROW>
            all_res(end+1)     = res_j;  %#ok<AGROW>
        end
    end
    
    info = struct();
    info.all_patches = all_patches;
    info.all_xi      = all_xi;
    info.all_eta     = all_eta;
    info.all_res     = all_res;
    info.ambiguous   = numel(all_patches) > 1;
    
    if isempty(all_patches)
        xi    = NaN;
        eta   = NaN;
        patch = NaN;
        warning('El punto (%.4f, %.4f) no pertenece a ning?n parche.', x, y);
    else
        xi    = all_xi(1);
        eta   = all_eta(1);
        patch = all_patches(1);

        if info.ambiguous
            fprintf('Aviso: el punto pertenece a %d parches: %s\n', ...
                    numel(all_patches), mat2str(all_patches));
        end
    end
    
end




function [xi, eta, conv, res] = invert_point(geom_j, x, y, xi0, eta0)
    % geometry_patch: la estructura geometry correspondiente a un solo parche
    % x, y: coordenadas del punto f?sico objetivo
    % xi0, eta0: estimaci?n inicial (por defecto (0.5, 0.5))

    if nargin < 4
        xi0 = 0.5;
        eta0 = 0.5;
    end

    tol = 1e-10;
    maxit = 50;
    conv= false;
    xi = xi0;
    eta = eta0;

% Obtener manejadores de funci?n (seg?n la estructura de GeoPDEs)
    F  = geom_j.map;
    DF = geom_j.map_der;
    for it = 1:maxit
        uv = [xi; eta];
        F_val = F(uv);
        J_val = DF(uv);
        r = [x; y] - F_val;
        res = norm(r);
        if res < tol
            conv=true;
            break;
        end
        delta = J_val \ r;
        xi  = xi  + delta(1);
        eta = eta + delta(2);
        % Limitar al dominio param?trico [0,1]^2 (opcional, evitar salirse)
        xi  = min(max(xi,  0), 1);
        eta = min(max(eta, 0), 1);
    end
    if it == maxit
        warning('Newton no convergi?, residual = %g', norm(r));
    end
end