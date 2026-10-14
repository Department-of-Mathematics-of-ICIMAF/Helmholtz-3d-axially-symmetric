function val = geometry_funct(Z, Z1, Z2, Z3)
    %% Funcion que da valores fijos segun la geometria (soportando escalares y matrices)
    
    % 1. Inicializamos 'val' del mismo tama?o que 'Z' lleno de ceros (caso por defecto)
    val = zeros(size(Z)); 
    
    % 2. Asignamos los valores usando mascaras logicas (indexacion vectorial)
    val(Z >= 0 & Z < Z1) = 1;
    val(Z >= Z1 & Z < Z2) = 2;
    val(Z >= Z2 & Z < Z3) = 3;
    val(Z >= Z3 ) = 3;
    
    % Nota: Los valores que no cumplen ninguna condicion ya se quedan en 0.

end