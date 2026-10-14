function[rho_dcoef,c_dcoef,k_dcoef,knots]=ShoenbergparaPennes_lentecurva(n,k,z1,z2,z3,z_min,r, show_figure)

% Esta funcion calcula un aproximante de Shoenberg para 3 funciones:
% la densidad rho(x,z), el calor especifico c_s(x,z) z la conductividad 
% termica k_s(x,z). Estas funciones son constantes para cada capa de tejido
% Por tanto son funciones constantes por pedazos para un tejido de 3 capas:
% piel, grasa z musculo 

% Entrada
% n,k- dimension z orden del espacio spline al que pertenece las funciones 
%    aproximantes de Shoenberg ( que solamente dependen de z)
% z1,z2,z3 - valores de z en (0,r) que definen las capas de tejido
% z_min- valor de z para el punto sobre la lente curvas con x=0
% r - radio del semicirculo ( 0< z1 < z2 < z3 < r)

% Salida

% rho_d (densidad),  c_d calor especifico,  k_d conductividad termica 
% En este problema estamos considerando que 
% para  0<= z <= 0.010     (zona de la piel)    rho_d=1200,   c_d=3590,  k_d=0.23
% para 0.010 < z <= 0.022  (zona de la grasa)   rho_d=950,    c_d=2670,  k_d=0.19
% para 0.022 < z <= 0.0355 (zona de musculo)    rho_d=1050,   c_d=3640,  k_d=0.55
% para 0.0355 < z <= r,    (zona del hueso)     rho_d=1975,   c_d=1250,  k_d=2.30



% El aproximante de Shoenberg para la densidad rho es: 
%  Vrho(z)=sum_{i=1}^n rho(tstar_i)B_i^k(z) 
gfun=@(z)density(z,r,z_min);
[rho_dcoef,knots]=ShoenbergSpline(n,k,gfun,z_min,r);
 
if show_figure ==1
    % % Evaluamos los B-splines en una secuencia fina de puntos en [a,b]
    tau= linspace(z_min, r, 10000);
    % Evaluamos el aproximante de Shoenberg
    Bmat = spcol(knots,k,tau);
    Vrho=Bmat*rho_dcoef;
    gtau=gfun(tau);
    figure(201), hold on
    plot(tau,gtau,'.-r')
    plot(tau,Vrho,'-b')
    title('Aproximante de Shoenberg para densidad \rho_d(0,z)'), axis tight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% El aproximante de Shoenberg para el calor especifico c_s es:
%  Vcs(z)=sum_{i=1}^n c_s(tstar_i)B_i^k(z) 
gfun=@(z)caloresp(z,r,z_min);
[c_dcoef,knots]=ShoenbergSpline(n,k,gfun,z_min,r);

if show_figure ==1 
    gtau=gfun(tau);
    figure(202), hold on
    plot(tau,gtau,'.-r')

    % Evaluamos el aproximante de Shoenberg
    %Bmat = spcol(knots,k,tau);
    Vcs=Bmat*c_dcoef;
    plot(tau,Vcs,'-b')
    title('Aproximante de Shoenberg para el calor especifico c_d(0,z)'), axis tight
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% El aproximante de Shoenberg para la conductividad termica k_s es:
%  Vks(z)=sum_{i=1}^n k_s(tstar_i)B_i^k(z) 
gfun=@(z)condter(z,r,z_min);
[k_dcoef,knots]=ShoenbergSpline(n,k,gfun,z_min,r);

if show_figure ==1
    gtau=gfun(tau);
    figure(203), hold on
    plot(tau,gtau,'.-r')

    % Evaluamos el aproximante de Shoenberg
    %Bmat = spcol(knots,k,tau);
    Vcs=Bmat*k_dcoef;
    plot(tau,Vcs,'-b')
    title('Aproximante de Shoenberg para la conductividad termica k_d(0,z)'), axis tight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function[rho]=density(z,r,z_min)
   nz=length(z);
   rho=zeros(nz,1);
   for i=1:nz
      if z_min <= z(i) && z(i) < 0
         rho(i)=1000; % water
      elseif 0 <= z(i) && z(i) <= z1 
         rho(i)=1200; % skin
      elseif z1 < z(i) && z(i) <= z2 
         rho(i)=950; % fat
      elseif z2 < z(i) && z(i) <= z3 
         rho(i)=1050; % muscle
      elseif z3 < z(i) && z(i) <= r         
%        rho(i)=1975; % bone
         rho(i)=1050; % muscle
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function[cs]=caloresp(z,r,z_min)
   nz=length(z);
   cs=zeros(nz,1);
   for i=1:nz
      if z_min <= z(i) && z(i) < 0
         cs(i)=4180; % water
      elseif 0 <= z(i) && z(i) <= z1 
         cs(i)=3590; % skin
      elseif z1 < z(i) && z(i) <= z2 
         cs(i)=2670; % fat
      elseif z2 < z(i) && z(i) <= z3 
         cs(i)=3640; % muscle
      elseif z3 < z(i) && z(i) <= r         
   % %         cs(i)=1250; % bone
         cs(i)=3640; % muscle
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function[ks]=condter(z,r,z_min)
   nz=length(z);
   ks=zeros(nz,1);
   for i=1:nz
      if z_min <= z(i) && z(i) < 0
         ks(i)=0.59; % water
      elseif 0 <= z(i) && z(i) <= z1 
         ks(i)=0.23; % skin
      elseif z1 < z(i) && z(i) <= z2 
         ks(i)=0.19; % fat
      elseif z2 < z(i) && z(i) <= z3 
         ks(i)=0.55; % muscle
      elseif z3 < z(i) && z(i) <= r         
% %         cs(i)=2.3; % bone
         ks(i)=0.55; % muscle
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
