function[b]=eval_q(space, msh, int_dofs, functionq)
% Esta fucion calcula el lado derecho
% 
% Entrada
% space - objeto de clase sp_scalar que representa el espacio de aproximacion spline 
% msh - objecto que contiene diferentes campos relativos a los puntos de ruptura, 
%      los nodos de la cuadratura y el mapeo que parametriza el dominio fisico
% int_dofs- vector de orden n0 que contiene los indices de las funciones
%            basicas que se anulan en la frontera ( o sea I0)
% Salida
% b - lado derecho (es de dimencion n0)
% 

functionb_rho = @(rho, z) rho.*functionq(rho,z);
b = op_f_v_mp(space, msh,functionb_rho); % Vector de tamanio Ndof
b = b(int_dofs); % Este es el vector b_s, s en I0
