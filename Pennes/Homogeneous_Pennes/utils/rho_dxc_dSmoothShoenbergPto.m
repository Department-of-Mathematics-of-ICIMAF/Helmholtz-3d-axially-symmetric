function[rho_dxc_drhoz]=rho_dxc_dSmoothShoenbergPto(rho,z,ksp,knots,rho_dcoef,c_dcoef)
% Esta funcion evalua en una matriz de puntos (rho,z) la Aproximacion Vcs(rho,z)de la 
% funcion c_s(rho,z) que describe el calor especifico en  un tejido multicapas. 
% Tambien evalua en (rho,z) la Aproximacion Vrho(rho,z)de la 
% funcion de densidad rho(rho,z). Estas funciones aparecen en la ecuacion de 
% Pennes z en la formulacion IgA forman parte del integrando para calcular 
% la matriz de masa
%  Entrada
% rho,z - matrices que definen los puntos donde se evalua Vks(rho(i,j),z(i,j))
% ksp, knots- orden z  sucesion de nodos de los splines de Shoenberg que se emplearon
%        para aproximar la conductividad termica kc(rho,z)
%  rho_dcoef- coeficientes del spline
%     Vrho_d(z)=sum_{i=1}^n rho_dcoef_i*B_i^k(z) - spline que aproxima a rho_d(z) (rho_d(rho,z) no depende de rho)

%  c_dcoef- coeficientes del spline
%     Vc_d(z)=sum_{i=1}^n c_dcoef_i*B_i^k(z) - spline que aproxima a c_d(z) (c_d(rho,z) no depende de rho)

%  Salida
% rho_dxc_drhoz - matriz del mismo tamanio que las matrices rho, z que contiene los
%       valores rho_dxc_drhoz(i,j)=Vrho_d((z(i,j)))*Vc_d(z(i,j))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,m]=size(rho);
rho_dxc_drhoz=zeros(n,m);
for i=1:n
  for j=1:m
    % Evaluamos los B-splines B_i^k(z) en los ptos z(i,:) de la fila i de z 
    Bmat = spcol(knots,ksp,z(i,j));
    % Evaluamos el spline Vcs en los ptos z(:,j) 
    % Vcs=Bmat*c_dcoef;
    c_drhoz=Bmat*c_dcoef;
    % Evaluamos el spline Vrho en los ptos z(:,j) 
    % Vrho=Bmat*rho_dcoef;
    rho_drhoz=Bmat*rho_dcoef;
    rho_dxc_drhoz(i,j)=rho_drhoz*c_drhoz;
  end
end
