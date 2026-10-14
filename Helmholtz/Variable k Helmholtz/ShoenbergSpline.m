function[gtstar,knots]=ShoenbergSpline(n,k,gfun,a,b)
% Esta funcion calcula la funcion Spline de Shoenberg Vg que aproxima a una
% funcion g conocida.
% Dada una sucesion de nodos t={t_1<=t_2<= ... <= t_{n+k}}
% el spline Vg se escribe como
%
% Vg(x)=sum_{i=1}^n g(tstar_i)B_i^k(x)
%
% donde B_i^k(x) son los B-splines de orden k para la sucesion de nodos t
% y tstar son puntos [t_1,t_{n+k}] que se calculan como promedio de los 
% nodos, o sea
% tstar_i=(t_{i+1}+...+t_{i+k-1})/(k-1)
% El spline Vg tiene buenas propiedades como aproximante:
% - si g es no negativa, entonces Vg tambien lo es
% - si g es monotona, entonces Vg tambien lo es
% - si g es convexa, entonces Vg tambien lo es
% - Vg aproxima localmente a g, pues el valor de Vg en [t_i,t_{i+1}]
% solamente depende de los valores de g en los puntos
% tstar_{i-k+1},...,tstart_{i}. Por ejemplo si los puntos 
% (tstar_{i-k+j},g(tstar_{i-k+j}))  j=1,...,k
% estan en una recta entonces Vg coincide con esa recta en [t_i,t_{i+k}]

% Entrada
% n - dimension del espacio spline
% k - orden del spline
% gfun- nombre de una funcion Matlab que permite evaluar la funcion g
% [a,b]- intervalo donde se construye el aproximante

% Salida

step=(b-a)/(n-(k-1));
% Puntos de ruptura
breakpts=a:step:b;
% Construccion de la sucesion de nodos 
knots= [a*ones(1,k-1) breakpts b*ones(1,k-1)];

% Calculo del promedio de ls nodos . Estos puntos se usan para definir la
% funcion spline Vg
tstar = aveknt(knots,k);

% Evaluacion de la funcion g en los ptos tstar
gtstar=gfun(tstar);

% Graficamos la funcion
%  Vg(x)=sum_{i=1}^n g(tstar_i)B_i^k(x)
%



