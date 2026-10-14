function[k_drhoz]=k_dSmoothShoenbergPto(rho,z,ksp,knots,k_dcoef)
% Esta funcion evalua en una matriz de puntos (rho,z) la Aproximacion Vk_d(rho,z)de la 
% funcion k_d(rho,z) que describe el coeficiente de conductividad termica de 
% un tejido multicapas. Este aparece en la ecuacion de Pennes y en la formulacion 
% IgA forma parte del integrando para calcular la matriz de rigidez
%  Entrada
% rho,z - matrices que definen los puntos donde se evalua Vk_d(rho(i,j),z(i,j))
% ksp, knots- orden y sucesion de nodos de los splines de Shoenberg que se emplearon
%        para aproximar la conductividad termica kc(rho,z)
%  k_dcoef- coeficientes del spline
%     Vk_d(z)=sum_{i=1}^n k_dcoef_i*B_i^k(z) - spline que aproxima a k_d(z) (k_d(rho,z) no depende de rho)
%  Salida
% k_drhoz - matriz del mismo tamanio que las matrices rho, z que contiene los
%       valores k_drhoz(i,j)=Vk_d(z(i,j))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,m]=size(rho);
k_drhoz=zeros(n,m);
for i=1:n
  for j=1:m
    % Evaluamos los B-splines B_i^k(z) en los ptos z(i,:) de la fila i de z 
    Bmat = spcol(knots,ksp,z(i,j));
    % Evaluamos el spline Vks en los ptos z(:,j) 
    % Vks=Bmat*ccoef;
    k_drhoz(i,j)=Bmat*k_dcoef;
  end
end