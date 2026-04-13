%mainHelmholtzMultipatch_3Daxial
% VHM marzo 2026
% programa para resolver la ecuacion de Helmholtz en un dominio 3D con 
% simetria axial, lente curva y k constante.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('------------------------------------------------------------------------\n')
fprintf('Solucion IgA de la ecuacion de Helmholtz 3D (simetria axial) en un problema\n')
fprintf('con lente curva \n')
fprintf('-------------------------------------------------------------------------\n')

% Datos del problema
r=0.133;
a=0.01;
c=1500;
fr=1e6;
%zc=input('zc=');
zc=r/2;
rl=sqrt(a^2+zc^2);
kval=2*pi*fr/c;
C=5.0e5;
fprintf('Wave number k=%e\n',kval);
fprintf('r=%e\n',r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcular la parametrizacion Multiparche del dominio con 3 parches 

alpha = 0.7698;
Q = [0.0217,0.0385];
% Parametrizacion Multiparche VHM con nodos simples
[geometry,boundaries,interfaces,boundary_interfaces]=ParametMultipatch_mediodominio(r,a,rl,alpha,Q,zc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construir el espacio spline al que pertenece la solucion 
% 
%n1=40;m1=24;n3=20;
%n1=80;m1=48;n3=40;
n1=100; m1=300;  n3=400; % N=308502 
%n1=50; m1=500;  n3=500; 
fprintf('-----------------------------------------------------------------\n')
fprintf('La dimension del espacio de aproximacion se calcula en terminos de\n')
fprintf('n1=%i,m1=%i,n3=%i\n',n1,m1,n3)

[space,msh]=approxspace_multipatch(geometry,boundaries,interfaces,boundary_interfaces,n1,m1,n3);
Ndof = space.ndof

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo de las matrices IgA

% Area exacta de Omega
% beta=atan(a/zc);
% Aexacta=pi*r^2/4+beta*rl^2/2-a*zc/2;
% fprintf('Area exacta de Omega=%5e\n', Aexacta);

u=radiationmatricesHelmholtz_3Daxial(space,msh,kval,C,r,alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Postprocesamiento grafico\n')
% Post procesamiento grafico
npoints=1000;
zp=linspace(0, 1, npoints);
pts = {0, zp};

% Evaluar los 2 parches de la solucion aproximada a lo largo del eje z
[up1, Fzp1] = sp_eval (u(space.gnum{1}), space.sp_patch{1}, geometry(1), pts);
Zp1=deal (squeeze(Fzp1(2,:,:)));
[up2, Fzp2] = sp_eval (u(space.gnum{2}), space.sp_patch{2}, geometry(2), pts);
Zp2=deal (squeeze(Fzp2(2,:,:)));

figure(400),grid on, hold on
plot(Zp1,real(up1),'LineWidth',2,'MarkerSize',6),grid on
plot(Zp2,real(up2),'LineWidth',2,'MarkerSize',6),grid on
title ('Parte Real de la solucion de Helmholtz sobre el eje z')

figure(500),grid on, hold on
plot(Zp1,imag(up1),'LineWidth',2,'MarkerSize',6),grid on
plot(Zp2,imag(up2),'LineWidth',2,'MarkerSize',6),grid on
title ('Parte imaginaria de la solucion de Helmholtz el eje z')

figure(600),grid on, hold on
plot(Zp1,abs(up1),'LineWidth',2,'MarkerSize',6),grid on
plot(Zp2,abs(up2),'LineWidth',2,'MarkerSize',6),grid on
title ('Valor absoluto de la solucion de Helmholtz sobre el eje z')

saveas(gcf,'multi.fig')
















