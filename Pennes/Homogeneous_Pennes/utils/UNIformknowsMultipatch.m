function[knots,zeta]=UNIformknowsMultipatch(geometry,degree,regularity,n1,m1,n3,npatch)
% Esta funcion calcula los nodos UNIFORMES de los 4 parches que definen el espacio de aproximacion
% para la geometria del semicirculo con ombligo

% Entrada
% geometry- estructura que describe el mapeo multipatch
% degree,regularity- grado y numero de derivadas continuas del espacio de funciones 
%                    cuyos nodos y ptos de ruptura se calculan en esta
%                    funcion y que se utilizara para aproximar la solucion
%                    de la pde.
% n1,m1 - Definen la dimension n1*m1 del espacio de funciones restringido a Omega^(1)
%         o sea el numero de nodos del espacio en Omega^(1) es n1+4 en la direccion xi y 
%        m1+4 en la direccion eta
% n3- Definen el numero de nodos en la direccion xi (n3+4) del espacio de funciones tanto
%      restringido a Omega^(3) como a Omega^(4)
%  El resto de las dimensiones quedan determinados en funcion de los anteriores, pues
%  n4=n3, m3=m4=m1
%  n2=n1, m2=n3
% npatch numero de parches
% Salida
% knots- celda 4X2- nodos de cada uno de los 4 espacios de aproximacion
%      Por ejemplo knots{1}{2} son los nodos teta del espacio restringido a Omega^(1)
%                  knots{3}{1} son los nodos txi del espacio restringido a Omega^(3)
% zeta - celda 4X2- puntos de ruptura de cada uno de los 4 espacios de aproximacion


nsub = [n1 m1; n1 n3; n3 m1; n3 m1];

for iptc = 1:npatch
% Para cada parche se insertan nodos mediante subdivision en cada direccion
% parametrica. De esta forma se construye el espacio de APROXIMACION LOCAL al que
% pertenece la solucion aproximada. Mas adelante se define el espacio de APROXIMACION
% GLOBAL a partir de los espacios de aproximacion locales de cada parche.
% Define the refined mesh, with tensor product structure
 [knots{iptc}, zeta{iptc}] = kntrefine (geometry(iptc).nurbs.knots, nsub(iptc,:)-4, degree, regularity);
end