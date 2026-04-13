function[kxy2]=ksmoothShoenbergPto(x,y,fr,ksp,knots,ccoef,alphacoef)
% Marzo 2022
% Esta funcion evalua en una matriz de puntos (x,y) la Aproximacion Vk(x,y)de la 
% funcion k(x,y) que aparece en la ecuacion de Helmholtz para el problema de 
% radiacion en un dominio semicircular de radio r. 
% ksmoothShoenberg es llamada por RadiationMatrixGeoPDEKcteaTramos cuando vamos a 
% calcular la matriz de masa y entra  como function handle de la funcion de 
% GeoDPE op_u_v_tp .
%
%  Entrada
% x,y - matrices que definen los puntos donde se evalua Vk(x(i,j),y(i,j))
% fr - frecuencia- se utiliza para definir el aproximante Vk(x,y)
% ksp, knots- orden y  sucesion de nodos de los splines de Shoenberg que se emplearon
%        para aproximar la velocidad c(x,y) del ultrasonido y la atenuacion
%        alpha(x,y)
%  ccoef- coeficientes del spline
%     Vc(y)=sum_{i=1}^n ccoef_i*B_i^k(y) - spline que aproxima a c(y) (c(x,y) no depende de x)
%  alphacoef- coeficientes del spline
%     Valpha(y)=sum_{i=1}^n alphacoef_i*B_i^k(y) - spline que aproxima a alpha(y) (alpha(x,y) no depende de x)
%  rhocoef-coeficientes del spline
%     Valpha(y)=sum_{i=1}^n rhocoef_i*B_i^k(y) - spline que aproxima a rho(y) (rho(x,y) no depende de x)
%  Salida
% kxy2 - matriz del mismo tamanio que las matrices x, y que contiene los
%       valores kxy(i,j)=(2*pi*fr/Vc(y(i,j))-i*Valpha(y(i,j)))^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cte=2*pi*fr;
[n,m]=size(x);
kxy2=zeros(n,m);
for i=1:n
  for j=1:m
   % Evaluamos los B-splines B_i^k(y) en los ptos y(i,:) de la fila i de y 
   Bmat = spcol(knots,ksp,y(i,j));
   % Evaluamos el spline Vc en los ptos y(:,j) 
   % Vc=Bmat*ccoef;
   cxy=Bmat*ccoef;
   if abs(cxy)<1e-6
       fprintf('cxy=%f\n',cxy)
   end
   % Evaluamos el spline Valpha en los ptos y(:,j) 
   % Valpha=Bmat*alphacoef
   alphaxy=Bmat*alphacoef;
   kxy2(i,j)=(cte/cxy-1i*alphaxy)^2;
%   kxy2(i,j)=(cte/cxy)^2;
  end
end
% fprintf('k2(x,y)\n')
% kxy2
% pause




