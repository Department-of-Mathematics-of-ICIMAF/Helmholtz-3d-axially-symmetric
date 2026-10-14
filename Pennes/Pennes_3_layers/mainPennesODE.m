% MainPennesODE
% Enero 2026
% 
% Programa para resolver la ecuacion de Pennes en un tejido de 3 capas.
%
% ??_d(x,y,z)c_d(x,y,z)(\partial T(x,y,z,t))/\partial t - k_d(x,y,z)\Delta T(x,y,z,t) = Q_e(x,y,z,t)
%
% (x,y,z) en Omega y t en [0,t_f]
% con Q_e(x,y,z,t)=Q(x,y,z)I_s(t)
% donde Q(x,y,z) es la intensidad del campo de presion acustica generado por el 
% transductor, la cual se calcula a partir del campo de presion acustica u(x,y,z) 
% que es la solucion del problema de Helmhotz 3D
% 
% Q(x,y,z) = ??_d(x,y,z)(|u(x,y,z)|^2)/(??_d(x,y,z)c_d(x,y,z))
% 
% donde:
%  ??_d - Dencidad del tejido
%  c_d - Calor especifico
%  k_d - Conductividad termica
%  ??_d - Atenuacion
% 
% Dado ?? > 0, I_s(t) es una funcion escalon C^2-continua, tal que I_s(t) = 1 
% para t ??? [0, t_s ??? 2??], I_s(t) = 0 para t ??? [t_s ??? ??, t_f] y I_s(t) es una 
% funcion decreciente para t ??? (t_s ??? 2??, t_s ??? ??).
% 
% Condicion inicial
% T (x, y, z, 0) = T_b
% donde T_b es constante en todo Omega.
% Condiciones de frontera
% T(x,y,z,t) = T_b para (x,y,z) en las fronteras "Abajo", "Derecha" y "Arriba" del 
% dominio Omega y para todo t en [0,tf]  (3)
% 
% Vamos a resolver el problema usando Analisis Isogeometrico en un dominio
% Omega que es una semi-esfera, que por simetria axial se reduce a un cuarto de circulo.
% Usando el m'etodo de Lineas el problema se reduce a un sistema ( de muchas!) 
% de ecuaciones diferenciales ordinarias ODE. Este sistema se escribe matricialmente  
% como
% M \alpha'(t) + K\alpha(t) = z(t)
%                          
% donde alpha(t)=(alpha_1(t),...,alpha_N(t)) es el vector de las funciones 
% incognitas en la representacion de la solucion aproximada Th en la base de B-splines, o sea
% Th(x,y,z,t) =sum_{i=1...N} alpha_i(t)* Phi_i(x,y,z)
% donde Phi_i(x,y,z) i=1...N es la base de B-splines
% El vector alpha_0(t) es de dimension n0 y tiene componentes alpha_r(t)
% con r en I_0, donde I_0 es el conjunto de indicies de las funciones
% basicas Phi_i(x,y,z) que se ANULAN para todo (x,y) en la frontera de Omega.

% M00 y K00 son matrices de orden n0, m'as precisamente bloques de las matrices 
% de masa M y de rigidez K.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
% Se agregan las carpetas a utilizar
addpath(genpath('./utils'));

theta = 3*pi/2;  % angulo de rotacion del solido
show_figure = 0; % 1 (true) or 0 (false) para mostrar los graficos de los resultados
% save_result = 0; % 1 (true) or 0 (false) para guardar los resultados en un txt

fprintf('------------------------------------------------------------------\n')
fprintf('Solucion IgA de la ecuacion de Pennes en un problema\n')
fprintf('donde se conoce la solucion exacta\n')
fprintf('------------------------------------------------------------------\n')

% Radio de la semiesfera
r=0.133; % radio de la semiesfera que define el dominio fisico
a=0.01 ; % semiapertura de la lente
rl = 0.0672; % Redio de curvatura de la lente
z_c= sqrt(rl^2 - a^2); % posicion del centro de la lente en el eje z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1=0.010; z2=0.032; z3=0.0675; % Grosor de las capas

zmin = z_c - rl; % Para la lente

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construir la parametrizacion F(chi,eta) del dominio fisico Omega
% F(chi,eta) se  puede parametrizar como una superficie reglada
[geometry,boundaries,interfaces,boundary_interfaces]=paramCuartCirculoLente(r, a, z_c, rl, show_figure);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construir el espacio spline al que pertenece la solucion 
fprintf('------------------------------------------------------------------\n')
% n1=4;m1=2;n3=2;
% n1=20;m1=12;n3=10;
% n1=40;m1=24;n3=20;
% n1 = 4;
n1 = 50;
m1=3*n1; n3=4*n1;
% n1=60;m1=48;n3=40;

[space,msh]=approxspace_multipatch(geometry,boundaries,interfaces,boundary_interfaces,n1,m1,n3);
Ndof = space.ndof;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculo de las matrices IgA
fprintf('----------------------------------------------------\n')
fprintf('Calculo de las matrices IgA basicas\n')
fprintf('orden de las matrices=%i\n',Ndof)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construir el aproximante de Shoenberg para las funciones constantes por
% capas que definen la densidad rho, el calor especifico c_s y la
% conductividad termica k_s
[rho_dcoef,c_dcoef,k_dcoef,knots]=ShoenbergparaPennes_lentecurva(50,3,z1,z2,z3,zmin,r,show_figure);

%
% Calcular las matrices IgA de masa y rigidez del espacio de aproximacion multipatch
fprintf('----------------------------------------------------------------\n')
fprintf('Calculo de las matrices IgA multipatch de masa y rigidez\n')
fprintf('----------------------------------------------------------------\n')
fprintf('tiempo para calcular las matrices\n')
tic
k_dfun=@(rho,z)k_dSmoothShoenbergPto(rho,z,3,knots,k_dcoef);
rho_dxc_dfun=@(rho,z)rho_dxc_dSmoothShoenbergPto(rho,z,3,knots,rho_dcoef,c_dcoef);


[M,K]=matricesIgAbasicas3D(space,msh,rho_dxc_dfun,k_dfun);

toc

%% Reducir las matrices M y K a los bloques correspondientes a los nodos que no estan
%% relacionados con la frontera de Dirichlet, que son los coeficientes que se tienen 
%% que calcular.
[drchlt_dofs, free_dofs] = select_drchlt_dofs(space);

M_0 = M(free_dofs, free_dofs);
K_0 = K(free_dofs, free_dofs);
K_d = K(free_dofs, drchlt_dofs);

%% Calculo de la condicion inicial
T_b = 35;
T_0 = @(rho,z) T_b;
% alpha_0_inicial = condinicial(space, msh, free_dofs, M_0,T_0);
alpha_0_inicial = T_b*ones(length(free_dofs),1);
n_d=length(drchlt_dofs);
alpha_d = T_b*ones(n_d,1);
%% Calculo del vector del lado derecho

% filenameQ=input('diga el nombre del fichero de valores de Q\n','s');
% filenameQ='Q.mat';
% filenameQ='QNdof189201_tejido3capas.mat';
filenameQ = 'Q_Dof189201_3capas_F1MHZ_C4e5.mat';
filename = sprintf('inputs//%s', filenameQ);
load(filename,'XYQ');
% lambda2 = 1.25*10^11;  % 1.25*10^11 factor para ajusatr temperatura alrededor de 50? en un tejido con 3 capas

%%% La escala lambda2 es el cociente entre
%%% la C que viene de Helmotz dividido por el nuevo C que se quiere para
%%% Pennes
% C_Helmotz = 1e6;
% C_Pennes = 4e5;
% lambda2 = (C_Pennes/C_Helmotz)^2;

lambda2 = 1;
Qval = scatteredInterpolant(XYQ(:,1),XYQ(:,2),lambda2*XYQ(:,3),'nearest');

q = eval_q(space, msh, free_dofs, Qval);  % no depende del tiempo

%% Calculo del sistema de ODEs
fprintf('-------------------------------------------------------------\n')
fprintf('Calculando la solucion del sistema ODE\n')
ATol = 1.0e-12; 
RTol = 1.0e-10; 
opts = odeset('RelTol',RTol,'AbsTol',ATol,'Mass',M_0,'Jacobian',-K_0);

% fprintf('El problema se resuelve en [0,tf]\n')
% tf=input('tf=');
tf = 48;
fprintf('El problema se resuelve en [0, %e]\n', tf)

tstop=30;
fprintf('La intensidad se detienen a los %i s\n',tstop)

t_in = 0;
N= 100;
step = (tf- t_in)/N;
% tspan = 0:step:tf;
tspan = [0 tf];

tic;
solver_name = 'ode23t';
[tiempo,alphaI0,S]= ode23t(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode23s';
% [tiempo,alphaI0,S]= ode23s(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode45';
% [tiempo,alphaI0,S]= ode45(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode15s';
% [tiempo,alphaI0,S]= ode15s(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),tspan, alpha_0_inicial, opts);

% solver_name = 'ode23';
% [tiempo,alphaI0,S]= ode23(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode113';
% [tiempo,alphaI0,S]= ode113(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);

% solver_name = 'ode23tb';
% [tiempo,alphaI0,S]= ode23tb(@(t, alpha_0) eval_lado_derecho(t, alpha_0, tstop, K_0, K_d, alpha_d, q),0:step:tf, alpha_0_inicial, opts);
solver_time = toc;

fprintf('Metodo de solucion %s\n',solver_name)
fprintf('-------------------------------------------------------------\n')
fprintf('Results of ODE solver\n')
fprintf('Number of successful steps = %i\n',S(1))
fprintf('Number of failed attemps = %i\n',S(2))
fprintf('Number of evaluations of the right hand side F(t,y)= %i\n',S(3))
fprintf('Number of evaluations of the Jacobian dF(t,y)/dy =%i\n',S(4))
fprintf('Number of LU decompositions = %i\n',S(5))
fprintf('Number of solutions of linear systems = %i\n',S(6))
fprintf('time to solve de ODE systems = %e\n',solver_time)
fprintf('-------------------------------------------------------------\n')

Nt=length(tiempo);

alpha(1:Nt, free_dofs) = alphaI0;
alpha(1:Nt, drchlt_dofs) = T_b*ones(Nt,n_d);



if show_figure == 0
    
    %% Graficos de la Solucion a lo largo del eje Z
    npoints=800;
    z_p = linspace(0,1, npoints);
    ptsejez = {0, z_p};
    max_value_g = 0;
% 
    t_5 = 5;
    [~,j_t_5]=min(abs(tiempo-t_5));
    [TempejezOmega1_tstop, Fzp]=sp_eval(alpha(j_t_5,space.gnum{1}),space.sp_patch{1},geometry(1),ptsejez);
    Z_tstop_Omega1=deal(squeeze(Fzp(2,:,:)));
    [TempejezOmega2_tstop, Fzp]=sp_eval(alpha(j_t_5,space.gnum{2}),space.sp_patch{2},geometry(2),ptsejez);
    Z_tstopOmega2=deal(squeeze(Fzp(2,:,:)));
    Z_tstop=[Z_tstop_Omega1' Z_tstopOmega2'];
    Tempejez_tstop=[TempejezOmega1_tstop TempejezOmega2_tstop];
    figure(500),grid on, hold on
    plot(Z_tstop,Tempejez_tstop,'-m','LineWidth',2,'MarkerSize',6),grid on
    max_value_g = max(max(Tempejez_tstop), max_value_g);
%     
    t_10 = 10;
    [~,j_t_10]=min(abs(tiempo-t_10));
    [TempejezOmega1_tstop, Fzp]=sp_eval(alpha(j_t_10,space.gnum{1}),space.sp_patch{1},geometry(1),ptsejez);
    Z_tstop_Omega1=deal(squeeze(Fzp(2,:,:)));
    [TempejezOmega2_tstop, Fzp]=sp_eval(alpha(j_t_10,space.gnum{2}),space.sp_patch{2},geometry(2),ptsejez);
    Z_tstopOmega2=deal(squeeze(Fzp(2,:,:)));
    Z_tstop=[Z_tstop_Omega1' Z_tstopOmega2'];
    Tempejez_tstop=[TempejezOmega1_tstop TempejezOmega2_tstop];
    plot(Z_tstop,Tempejez_tstop,'-k','LineWidth',2,'MarkerSize',6),grid on
    max_value_g = max(max(Tempejez_tstop), max_value_g);
% 
    t_15 = 15;
    [~,j_t_15]=min(abs(tiempo-t_15));
    [TempejezOmega1_tstop, Fzp]=sp_eval(alpha(j_t_15,space.gnum{1}),space.sp_patch{1},geometry(1),ptsejez);
    Z_tstop_Omega1=deal(squeeze(Fzp(2,:,:)));
    [TempejezOmega2_tstop, Fzp]=sp_eval(alpha(j_t_15,space.gnum{2}),space.sp_patch{2},geometry(2),ptsejez);
    Z_tstopOmega2=deal(squeeze(Fzp(2,:,:)));
    Z_tstop=[Z_tstop_Omega1' Z_tstopOmega2'];
    Tempejez_tstop=[TempejezOmega1_tstop TempejezOmega2_tstop];
    plot(Z_tstop,Tempejez_tstop,'-r','LineWidth',2,'MarkerSize',6),grid on
    max_value_g = max(max(Tempejez_tstop), max_value_g);
% 
    t_25 = 25;
    [~,j_t_25]=min(abs(tiempo-t_25));
    [TempejezOmega1_tstop, Fzp]=sp_eval(alpha(j_t_25,space.gnum{1}),space.sp_patch{1},geometry(1),ptsejez);
    Z_tstop_Omega1=deal(squeeze(Fzp(2,:,:)));
    [TempejezOmega2_tstop, Fzp]=sp_eval(alpha(j_t_25,space.gnum{2}),space.sp_patch{2},geometry(2),ptsejez);
    Z_tstopOmega2=deal(squeeze(Fzp(2,:,:)));
    Z_tstop=[Z_tstop_Omega1' Z_tstopOmega2'];
    Tempejez_tstop=[TempejezOmega1_tstop TempejezOmega2_tstop];
    plot(Z_tstop,Tempejez_tstop,'-g','LineWidth',2,'MarkerSize',6),grid on
    max_value_g = max(max(Tempejez_tstop), max_value_g);
% 
     t_29 = 30;
    [~,j_tstop]=min(abs(tiempo-t_29));
    [TempejezOmega1_tstop, Fzp]=sp_eval(alpha(j_tstop,space.gnum{1}),space.sp_patch{1},geometry(1),ptsejez);
    Z_tstop_Omega1=deal(squeeze(Fzp(2,:,:)));
    [TempejezOmega2_tstop, Fzp]=sp_eval(alpha(j_tstop,space.gnum{2}),space.sp_patch{2},geometry(2),ptsejez);
    Z_tstopOmega2=deal(squeeze(Fzp(2,:,:)));
    Z_tstop=[Z_tstop_Omega1' Z_tstopOmega2'];
    Tempejez_tstop=[TempejezOmega1_tstop TempejezOmega2_tstop];
    plot(Z_tstop,Tempejez_tstop,'-b','LineWidth',2,'MarkerSize',6),grid on
    max_value_g = max(max(Tempejez_tstop), max_value_g)
% 
    plot([z1,z1],[35,max_value_g],'--c','LineWidth',2)
    plot([z2,z2],[35,max_value_g],'--c','LineWidth',2)
% 
    legend('t=5 sec', 't=10 sec', 't=15 sec', 't=25 sec', 't=30 sec')
    xlabel('Temp.'); ylabel('Z');
    ax = gca;
    ax.FontSize = 12;
    ax.FontWeight = 'bold';
    
    hold off
    
    %% Graficos de la Solucion 3D
    npoints=100;
    pts = {linspace(0, 1, npoints), linspace(0, 1, npoints)};
    
    n_patch = space.npatch;
    
    min_value_g = 100;
    max_value_g = 0;
    
    geometry_function_hande = @(rho, z) geometry_funct(z, z1, z2, z3);
    % 1. Definir los colores exactos en formato RGB (valores de 0 a 1)
    azul      = [0.0, 0.4, 0.8]; % Lente
    carmelita = [0.7, 0.5, 0.3]; % Piel (marr?n/carmelita claro)
    amarillo  = [1.0, 0.9, 0.2]; % Grasa
    rojo      = [0.9, 0.3, 0.5]; % M?sculo

    % 2. Crear la matriz de la paleta (el orden debe coincidir con tus Ticks 0, 1, 2, 3)
    mi_paleta = [
        azul;        % Posici?n 0
        carmelita;   % Posici?n 1
        amarillo;    % Posici?n 2
        rojo         % Posici?n 3
    ];
    for patch = 1:n_patch    
% %         j = 60;
        t_eval = 30;
        [~,j] = min(abs(tiempo-t_eval));
        [Taprox, F] = sp_eval (alpha(j,space.gnum{patch}), space.sp_patch{patch}, geometry(patch), pts);
        [Rho, Z]  = deal(squeeze(F(1,:,:)), squeeze(F(2,:,:)));
%     % % % % %  Figura con la solucion calculada en un tiempo fijo % % % % %
        min_value = min(Taprox);
        if length(min_value)>1
            min_value = min(min_value);
        end
        max_value = max(Taprox);
        if length(max_value)>1
            max_value = max(max_value);
        end
        min_value_g = min(min_value, min_value_g);
        max_value_g = max(max_value, max_value_g);    
%         figure(600), hold on, grid off, view(3)
%         surf(Rho, Z,  Taprox, 'EdgeColor', 'interp')
%         title(['Approximated solution of Pennes equation in t= ', num2str(tiempo(j)), 'sec']) 
%         zlabel('U')
%         xlabel('\rho')
%         ylabel('Z')
%         colorbar('Ticks',linspace(min_value_g,max_value_g,6));
%         hold off

        figure(602), hold on,grid on
        plot_3D_cylinder_1(Rho, Z, Taprox, theta, min_value_g, max_value_g, ['Temperatura para t = ', num2str(tiempo(j)), 'sec'],z1,z2,z3);
        title(['Superficie rotada para el timepo t=', num2str(tiempo(j)), 'sec']);
        %%% Grafico las curvas de los bordes
        hold on
        x_rho = r*sin(0:0.001:pi/2);
        x_z = r*cos(0:0.001:pi/2);
        x_zero = zeros(1,length(x_rho));
        plot3(x_rho,x_zero,x_z, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'w'); % 
        plot3(x_zero, -x_rho,x_z, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'w'); % 
        hold off
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
        hold off
        
        
        figure(603), hold on,grid on
        surf(Rho, Z, Taprox,'EdgeColor', 'flat')
        surf(-Rho, Z, Taprox,'EdgeColor', 'flat')
        title(['Superficie rotada para el timepo t=', num2str(tiempo(j)), 'sec']);
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
        cb = colorbar;
        cb.Limits = [min_value_g, max_value_g];
        xlabel('X'); ylabel('Z');
        hold off
        
%         figure(800), hold on
%         ZZ = Qval(Rho,Z);
%         max_Q = max(max(ZZ))
%         surf(Rho, Z, ZZ, 'EdgeColor', 'interp')
%         colorbar
%         title('Valor de la Intencidad Q de entrada.');
%         hold off
         
%         figure(801), hold on%, grid on, view(2)
%         surf(Rho, Z, geometry_function_hande(Rho, Z),'EdgeColor', 'flat')
%         plot_3D_cylinder(Rho, Z, geometry_function_hande(Rho, Z), theta);
%         colorbar
%         title (['Exact solution of Pennes equation in t= ', num2str(tiempo(j)), 'sec']),
%         zlabel('U')
%         ylabel('Z')
%         xlabel('\rho')
%         colormap(mi_paleta);
%         cb = colorbar;
%         cb.TickLabels = {'lent','skin','fat','muscle'};
%         cb.Ticks = [0.25, 1, 1.8, 2.6];
%         grid on
%         hold off

    end
end