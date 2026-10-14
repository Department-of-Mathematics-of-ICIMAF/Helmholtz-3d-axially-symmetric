function[alpha0]=condinicial(space, msh, int_dofs, M_0,T_0)
% Esta funcion calcula los coeficientes de la condicion inicial
% Entrada
% space - objeto de clase sp_scalar que representa el espacio de aproximacion spline 
% msh - objecto que contiene diferentes campos relativos a los puntos de ruptura, 
%      los nodos de la cuadratura y el mapeo que parametriza el dominio fisico
% int_dofs- vector de orden n0 que contiene los indices de las funciones
%            basicas que se anulan en la frontera ( o sea I0)
% 
% Salida
% alpha0 - coeficientes de la aproximacion de las condiciones iniciales (es de dimencion Int_dof)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
functionT0_rho = @(rho, z) rho.*T_0(rho,z);
b=op_f_v_mp(space, msh,functionT0_rho); % Vector de tamanio Ndof
b = b(int_dofs);    % Este b es de tamanio Int_dof

% Resolver el sistema de ecuaciones lineales
%   M_0*alpha0=b_0
% para obtener los coeficientes alpha0 de la funcion que aproxima 
% a T (�?, z, 0)

alpha0 =M_0\b; % Este es el vector alpha_s, s en I0
