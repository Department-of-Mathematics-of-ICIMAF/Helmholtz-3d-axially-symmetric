function[space,msh]=approxspace_multipatch(geometry,boundaries,interfaces,boundary_interfaces,n1,m1,n3)
% Esta funcion construye el espacui de aproximacion basada en funciones
% B-spline multipatch al que pertenece la solucion de la ecuacion de
% Pennes
%
% Entrada
% geometry- estructura que describe el mapeo multipatch, parametro a la salida de  mp_geo_load
% boundaries,interfaces,boundary_interfaces- parametros de salida de mp_geo_load
% n1,m1,n3- Parametros que definen la dimension del espacio de aproximacion.
%           El espacio esta compuesto por 4 subespacios V^j, cada uno definido
%           sobre un subdominio Omega^j, j=1,...,4
%           La dimension de V^j es nj*mj donde n1, m1 y n3 se escogen libremente mientras que
%           el resto de las dimensiones quedan determinados en funcion de los anteriores, pues
%           n4=n3, m3=m4=m1
%           n2=n1, m2=n3
fprintf('----------------------------------------------------------------\n')
fprintf('Contruyendo el espacio de aproximacion para cada parche\n')

npatch = numel (geometry);
% Para cada parche tenemos un vector de 2 componentes, la primera es el numero
%de nodos que se insertan entre cada par de puntos de ruptura consecutivos
%en la direccion horizontal y la segunda en la direccion vertical

nquad = [4 4];
regularity = [2 2];
degree = [3 3];
msh = cell (1, npatch); 
sp = cell (1, npatch);


% Construccion de las sucesiones UNIFORMES de nodos y puntos de ruptura de cada parche
 [knots,zeta]=UNIformknotsMultipatch(geometry,degree,regularity,n1,m1,n3,npatch);

for iptc = 1:npatch
%
% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  local_msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));

% sp: object representing the discrete function space
   sp{iptc} = sp_bspline (knots{iptc}, degree, local_msh{iptc});
end

% Aqui se construye el mesh msh y el espacio spline globales, que incluyen 
% la informacion de los 4 parches. 
msh = msh_multipatch (local_msh, boundaries);
space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
% space.sp_patch{1}.space_type nos dice si el tipo de espacio del parche 1
% es 'spline' o 'NURBS'. En este problema los 4 parches son 'spline'

% space.sp_patch{1}.ndof - nos dice la dimension del espacio Vh
% correspondiente al parche 1


