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
fprintf('tiempo para construir el espacio de aproximacion para cada parche\n')

npatch = numel (geometry);
% Para cada parche tenemos un vector de 2 componentes, la primera es el numero 
% de subdivisiones horizontales y la segunda es el numero de subdivisiones verticales
% subdivisiones horizontales. Observar que el numero de subdivisiones en las fronteras 
% comunes a dos parches debe ser el mismo. 
%nsub = [10 6; 10 5; 5 6; 5 6];
%nsub = [20 12; 20 10; 10 12; 10 12];
%nsub = [40 24; 40 20; 20 24; 20 24]; % fue el que use para resolver el problema del paper!
%nsub = [80 48; 80 40; 40 48; 40 48];
%nsub = [100 60; 100 50; 50 60; 50 60];
%nsub = [200 120; 200 100; 100 120; 100 120];
nquad = [4 4];
regularity = [2 2];
degree = [3 3];
msh = cell (1, npatch); 
sp = cell (1, npatch);
tic

% Construccion de las sucesiones UNIFORMES de nodos y puntos de ruptura de cada parche
 [knots,zeta]=UNIformknowsMultipatch(geometry,degree,regularity,n1,m1,n3,npatch);

for iptc = 1:npatch
% Esto YA NO se usa!
% Para cada parche se insertan nodos mediante subdivision en cada direccion
% parametrica. De esta forma se construye el espacio de APROXIMACION LOCAL al que
% pertenece la solucion aproximada. Mas adelante se define el espacio de APROXIMACION
% GLOBAL a partir de los espacios de aproximacion locales de cada parche.
% Define the refined mesh, with tensor product structure
%  [knots{iptc}, zeta{iptc}] = ...
%         kntrefine (geometry(iptc).nurbs.knots, nsub(iptc,:)-1, degree, regularity);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

toc
