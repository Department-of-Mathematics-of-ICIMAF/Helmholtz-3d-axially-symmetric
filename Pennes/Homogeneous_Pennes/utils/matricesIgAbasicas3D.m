function[mass_mat,stiff_mat]=matricesIgAbasicas3D(space,msh,M_funct,K_funct)
% Enero 2026, VHM
% Esta funcion calcula las matrices M, K, de la discretizacion IgA de 
% la ecuacion de Pennes 3D que se proyecta a 2D usando simetria axial
% El dominio Omega es 1/4 de circulo.
%
% Entrada
% space - objeto de clase sp_scalar que representa el espacio de aproximacion spline 
% msh - objecto que contiene diferentes campos relativos a los puntos de ruptura, 
%      los nodos de la cuadratura y el mapeo que parametriza el dominio fisico
% Nof=n*m numero total de grados de libertad
% M_funct - funciones que se multiplican dentro de la matriz de masa
% K_funct - funciones que se multiplican dentro de la matriz de rigidez
%
% Salida
% K- matriz de rigidez de orden Ndof cuyo elemento (i,j) es 
%    k_{i,j}=2*pi int int rho*(dBi/drho,dBi/dz)*(dBj/drho,dBj/dz)' drho dz
% M - matriz de masa de orden Ndof cuyo elemento (i,j) es
%    m_{i,j}=2*pi int int rho* Bi* Bj drho dz
% int_dofs-indices I0 de las funciones basicas que se anulan sobre la frontera de Dirichlet    
% drchlt_dofs-indices Ig de las funciones basicas que NO se anulan sobre la frontera

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculo de la matriz de rigidez de orden Ndof x Ndof. Se almacena de forma sparse
frho = @(rho,z) rho;

k_f= @(rho,z) K_funct(rho,z).*frho(rho,z);
stiff_mat= op_gradu_gradv_mp (space, space, msh, k_f);

% Calculo de la matriz de masa de orden Ndof x Ndof. Se almacena de forma sparse
m_f= @(rho,z) M_funct(rho,z).*frho(rho,z);
mass_mat = op_u_v_mp (space,space, msh, m_f);
                                                 