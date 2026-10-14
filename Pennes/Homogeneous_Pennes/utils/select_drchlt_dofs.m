function [I_Dir, I_int] = select_drchlt_dofs(space)
  % Numeracion local de las fronteras de un parche
  % frontera 1 - Oeste
  % frontera 2 - Este
  % frontera 3 - Sur
  % frontera 4 - Norte
  %
  % Para el problema 3D con simetria axial y dominio definido por 1/4 de
  % circulo union el otro cuarto de circulo que define la lente curva, tenemos 
  % definida una condicion de Dirichlet en Robin homogenea en 4 subfronteras 
  % Como se numeran las fronteras??
  % Respuesta: Siguiendo el orden de los parches. En este problema tenemos 3
  % parches P1 en el centro abajo y a la izquierda de Omega, P2 sobre P1 y P3 
  % a la derecha de P1 y P2. Entonces tenemos 6 segmentos de frontera: la
  % primera es la oeste de P1, la segunda es la Sur de P1, despues la Oeste de 
  % de P2 y la Norte de P2 en ese orden y finalmente la Este y la Sur de P3 en 
  % ese orden.

  % Indices globales de las funciones basicas que no se anulan sobre cada una
  % de las 6 subfronteras en que se particiona la frontera de Omega
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % La solucion aproximada  % T(rho,z,t) de la ecuacion de Pennes 3D con simitria axial 
  % se escribe como
  %  T(rho,z,t)= sum_{j=1}^Ndof alpha_j(t)*Phi_j(rho,z)
  % donde Phi_j(rho,z) son los B-splines producto tensorial y alpha_j(t) son 
  % funciones incognitas que se calculan exigiendo que se satisfaga la
  % formulacion de Galerkin. 
  %
  % Los "coeficientes" alpha_j(t), j=1,...,N  se dividen en 2 conjuntos: 
  % {1,...,N}= (I_Dir) union (I_noDir), donde I_Dir es el conjunto de los indices 
  % (globales) de las funciones que NO se anulan sobre la frontera donde se ha definido
  % una condicion de Dirichlet ( I_noDir es el resto de los indices). A su
  % vez I_Dir se puede escribir como 
  %   I_Dir=(IDir_g) union (IDir_0)
  % donde IDir_g es el conjunto de los indices globales de las funciones
  % basicas que no se anulan sobre la frontera de la lente, o sea IDir_g=iglobalboun2.
  % IDir_0 es el conjunto de los indices globales de las funciones basicas que 
  % no se anulan sobre el resto de la frontera de Dirichlet (donde la condicion
  % de Dirichlet es homogenea), o sea 
  % IDir_0=(iglobalboun4) union (iglobalboun5) union (iglobalboun6)

  drchlt_global_sides = [2 4 5 6];
  I_Dir = [];
  I_int = [];

  % Ciclo sobre las fronteras de Dirichlet para determinar los indices drchlt_dofs 
  % de las funciones que NO se anulan sobre la frontera de Dirichlet. Los
  % coeficientes con indices drchlt_dofs en la representacion de la solucion 
  % del problema en la base de B-splines se utilizan para garantizar que se 
  % cumpla la condicion de Dirichlet.
  for j= 1:length(drchlt_global_sides)
    I_Dir = union(I_Dir, space.boundary.dofs(space.boundary.gnum{drchlt_global_sides(j)}));
  end

  Ndof=space.ndof;
  I=1:Ndof;
  I_int=setdiff(I,I_Dir);
  
end