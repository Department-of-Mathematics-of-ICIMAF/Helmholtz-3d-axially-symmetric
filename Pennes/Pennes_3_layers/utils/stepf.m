function[sp]=stepf
% Esta funcion calcula el spline cubico que interpola valores de una
% funcion escalon.

stepdata=importdata('funcion step.txt');
%  figure(111), hold on
%  plot(stepdata(:,1),stepdata(:,2),'ob')

% funcion tipo escalon
%y=0.5*(1+tanh(-750*x))  x en [-0.0075,0.0075]
% Esta funcion no aproxima muy bien cerca de -0.0075 y 0.0075
% Por eso no la usamos finalmente
% x=stepdata(:,1);
% y=0.5*(1+tanh(-750*x));
% plot(x,y,'-b')

% Usamos en su lugar un spline cubico de interpolacion con nodos escogidos 
% como el promedio de los datos
sp=spapi(4,stepdata(:,1),stepdata(:,2));

% fnplt(sp,2,'r')
% fprintf('---\n')
