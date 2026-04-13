function[u]=IgALenteCurvo_Multipatch_kvariable(geometry,space,msh,CDir,kfun2,kfunR,knots,ccoef,alphacoef,rhocoef)
% Enero 2025
%
% Esta funcion utiliza el analisis isogeometrico para resolver la ecuacion
% de Helmholtz
% -du^2/dx^2 - du^2/dy^2 -k(x,y)^2*u =0
% en Omega= semicircle of radio r and center in (0,0) + transductor curvo
% con condiciones mixtas en la frontera de Omega:
%
% Robin
% du/dn +i*k*u=gr  en el semicirculo
%
% Neuman
% du/dn=0 en los segmentos de la base del semicirculo, complemento de la lente curva
%
% Dirichlet
% u=gD sobre la lente curva
%
% Entrada
%  geometry- estructura que describe el mapeo multipatch, parametro a la salida de  mp_geo_load
%  space -  multipatch space, formed by several tensor product spaces plus the connectivity (see sp_multipatch)
%  msh   -  multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%  CDir  -  constante de Dirichlet, amplitud del pulso
%  kfun2- nombre de la funcion Matlab para evualar la funcion constante por pedazos k(x,y)^2 
%        Esta se utiliza para calcular la matriz de Masa
%  kfunR- nombre de la funcion Matlab para evualar la funcion constante por pedazos k(x,y). 
%        Esta se utiliza para calcular la matriz de Robin
%
% Salida
% u - vector solucion del sistema de ecuaciones lineales, o sea los
%     coeficientes en la base de B-splines de la solucion del problema de
%     radiacion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensamblar la matriz A y calcular el lado derecho  b 
%fprintf('time for computing matrix and right hand side \n')
%fprintf('----------------------------------------------------------------\n')
tic

% Calcular las matrices de masa, rigidez y Robin
c_diff  = @(x, y) ones(size(x));

% Calculo de la matriz de rigidez
stiff_mat = op_gradu_gradv_mp (space, space, msh, c_diff);

% Calculo de la matriz de masa. OJO como k(x,y) NO es constante entonces
% el elemento M_{i,j} de la matriz de masa es 
% M_{i,j}= int int_ Omega (k(x,y)^2* Phi_i(x,y)*Phi_j(x,y) dx dy
% Para calcular esta matriz hay que entrarle a la funcion op_u_v_mp un 4to
% parametro que es una funcion Matlab que evalue a k(x,y)^2

% Calculo de la matriz de masa
mass_mat = op_u_v_mp(space,space, msh, kfun2);

% Inicialmente la matriz global mat es stiff_mat - mass_mat, despues
% hay que imponer las condiciones de frontera de Dirichlet y Robin
mat=stiff_mat - mass_mat;

Ndof = size(mat,1);
%--------------------------------------------------------------------------
%
% Calculo de la matriz de Robin.
%
% Numeracion de las fronteras de un parche
% frontera 1 - Oeste
% frontera 2 - Este
% frontera 3 - Sur
% frontera 4 - Norte
%
% Tenemos definida la condicion de Robin homogenea en  3 subfronteras 
% Como se numeran las fronteras??
% Respuesta: Siguiendo el orden de los parches. En este problema tenemos 4
% parches P1 en el centro abajo de Omega, P2 sobre P1, P3 a la izquierda de
% P1 y P4 a la derecha de P1. Entonces tenemos 6 segmentos de frontera: la
% primera es la Sur de P1, despues la Norte de P2, despues la Oeste y la
% sur de P3 (en ese orden) y finalmente la Este y la Sur de P4 (en ese
% orden)

% La matriz robin1_mat corresponde a la frontera global 3, que es la frontera 1 (Oeste) del parche 3 
 [rows_robin_mat1, cols_robin_mat1, val_robin_mat1] = op_u_v_tp(space.sp_patch{3}.boundary(1), ...
                                                      space.sp_patch{3}.boundary(1), msh.msh_patch{3}.boundary(1),kfunR);
% La matriz robin2_mat corresponde a la frontera global 2, que es la frontera 4 (Norte) del parche 2
 [rows_robin_mat2, cols_robin_mat2, val_robin_mat2] = op_u_v_tp(space.sp_patch{2}.boundary(4), ...
                                                      space.sp_patch{2}.boundary(4), msh.msh_patch{2}.boundary(4),kfunR);
% La matriz robin3_mat corresponde a la frontera global 5, que es la frontera 2 (Este) del parche 4
 [rows_robin_mat3, cols_robin_mat3, val_robin_mat3] = op_u_v_tp(space.sp_patch{4}.boundary(2), ...
                                                      space.sp_patch{4}.boundary(2), msh.msh_patch{4}.boundary(2),kfunR);
 NRobin1=space.sp_patch{3}.boundary(1).ndof; %numero de funciones basicas que no se anulan en la frontera de Robin 1 y por
                                          % tanto dimension de la matriz de Robin robin1_mat
 NRobin2=space.sp_patch{2}.boundary(4).ndof; %numero de funciones basicas que no se anulan en la frontera de Robin 2 y por
                                          % tanto dimension de la matriz de Robin robin2_mat
 NRobin3=space.sp_patch{4}.boundary(2).ndof; %numero de funciones basicas que no se anulan en la frontera de Robin 3 y por
                                          % tanto dimension de la matriz de Robin robin3_mat  
                                          
 % Ensamblando las matrices                                        
 robin1_mat=sparse(rows_robin_mat1, cols_robin_mat1, val_robin_mat1,NRobin1,NRobin1);
 robin2_mat=sparse(rows_robin_mat2, cols_robin_mat2, val_robin_mat2,NRobin2,NRobin2);
 robin3_mat=sparse(rows_robin_mat3, cols_robin_mat3, val_robin_mat3,NRobin3,NRobin3);

% Indices globales de las funciones basicas que no se anulan sobre cada una
% de las 6 subfronteras en que se particiona la frontera de Omega
iglobalboun1=space.boundary.dofs(space.boundary.gnum{1}); % frontera de abajo del parche 1, ombligo de Omega
iglobalboun2=space.boundary.dofs(space.boundary.gnum{2}); % frontera de arriba del parche 2, arco central del semicirculo 
iglobalboun3=space.boundary.dofs(space.boundary.gnum{3}); % frontera izquierda del parche 3, arco izquierdo del semicirculo
iglobalboun4=space.boundary.dofs(space.boundary.gnum{4}); % frontera de abajo del parche 3, segmento de recta
iglobalboun5=space.boundary.dofs(space.boundary.gnum{5}); % frontera de derecha del parche 4, arco derecho del semicirculo
iglobalboun6=space.boundary.dofs(space.boundary.gnum{6}); % frontera de abajo del parche 4, segmento de recta

% Finalmente actualizamos las filas y columnas de mat cuyos indices corresponden a los B-splines que no se 
% anulan sobre la frontera de Robin 
% mat = stiff_mat - mass_mat + i* robin_mat
mat(iglobalboun3,iglobalboun3)=mat(iglobalboun3,iglobalboun3)+1i*robin1_mat; % 
mat(iglobalboun2,iglobalboun2)=mat(iglobalboun2,iglobalboun2)+1i*robin2_mat; % 
mat(iglobalboun5,iglobalboun5)=mat(iglobalboun5,iglobalboun5)+1i*robin3_mat; % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Condicion de Dirichlet
dofs_dirich=iglobalboun1;
alpha1=CDir*ones(length(dofs_dirich),1);
u(dofs_dirich)=alpha1;

%------------------lado derecho---------

% Actualizacion del lado derecho despues de fijar las incognitas relacionadas con 
% la condicion de Dirichlet

% Grados de libertad interiores, o sea indices de las funciones B-spline
% que se anulan sobre la frontera de Dirichlet
int_dofs = setdiff (1:space.ndof, dofs_dirich);
rhs=zeros(Ndof,1);

rhs(int_dofs) = rhs(int_dofs) - mat(int_dofs, dofs_dirich) * u(dofs_dirich)';
fprintf('tiempo para calcular la matriz y el lado derecho\n')

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Aqui comienza el Algebra lineal
% Calcular la solucion del sistema lineal despues de eliminar las
% incognitas que se prefijaron para garantizar que se cumpla la condicion
% de Dirichlet
%fprintf('Time for solving the linear system\n')
fprintf('----------------------------------------------------------------\n')

fprintf('Seleccionar el metodo de solucion del sistema lineal Ax=b\n')
fprintf('Opcion 0 = metodo directo \n')
fprintf('Opcion 1 = GMRES con Complex Shiftex Laplacian Preconditioner\n')
var=input('Option=');
var=0;
fprintf('Resolviendo el sistema lineal ...\n')
tic
if var==0       
 u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);
elseif var==1
 u(int_dofs)=gmres_iluCSL(mat(int_dofs, int_dofs),rhs(int_dofs),mass_mat(int_dofs, int_dofs),kval);
else
 error('error selecting the option for solving linear system') 
end
fprintf('Tiempo para resolver el sistema lineal\n')
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XYQ]=evaluaQ(knots,ccoef,alphacoef,rhocoef,u,space,geometry);
% Experimentos variando la constante de Dirichlet CDir
filename=input('diga el nombre del fichero de valores de Q\n','s');

%filename='valoresdeQKvar_n1100m1300n3400_CDir0p25ycrover2.mat';
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p5ycrover2.mat';
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p75ycrover2.mat';
%filename='valoresdeQKvar_n1100m1300n3400_CDir1p0ycrover2.mat';

%filename='valoresdeQKvar_n1100m1300n3400_CDir1p0ycrover8.mat';

% Experimentos variando la posicion del centro de la lente y por tanto de
% su curvatura
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p5por10ala6ycrover2.mat';
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p5por10ala6ycrover4.mat';
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p5por10ala6ycrover6.mat';
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p5por10ala6LentePLANA.mat';

% Experimentos variando la semiapertura
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p5por10ala6semiaper5mm.mat';
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p5por10ala6semiaper10mm.mat';
%filename='valoresdeQKvar_n1100m1300n3400_CDir0p5por10ala6semiaper15mm.mat';

save(filename,'XYQ')



