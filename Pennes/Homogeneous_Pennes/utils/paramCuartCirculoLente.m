function[geometry,boundaries,interfaces,boundary_interfaces]=paramCuartCirculoLente(r, a, z_c, rl, show_figure)

% Esta funcion construye una parametrizacion F(xi,eta) del dominio del problema 
% 3D con lente CURVA cuando se utiliza la simetria axial. Este dominio es 2D y esta 
% definido por 1/4 de circulo de radio r al cual se le agrega el sector de circulo
% que define la lente. 
% 
% Entrada
% r- radio del circulo
% a- semiapertura de la lente
% z_c- posicion del centro de la lente en el eje z 
%
% Salida
% geometry- Estructura con diferentes campos creada por mp_geo_load
% boundaries - ???
% interfaces - ???
% boundary_interfaces - ???

alpha = 0.7698;
P_int=[0.0217,0.0385];
r_P_int=norm(P_int);
beta=acos(P_int(1)/r_P_int);

line_right1 = abline([a, 0],P_int);
line_right2 = abline(P_int,[r*cos(alpha), r*sin(alpha)]);
line_right3 = abline([a, 0],[r, 0]);

line_left1 = abline([0, -(rl-z_c)],[0, r_P_int]);
line_left2 = abline([0, r_P_int], [0, r]);

arc_top = nrbcirc(r,[0 0],alpha,pi/2);
aux=arc_top.coefs(:,1);
arc_top.coefs(:,1)=arc_top.coefs(:,3);
arc_top.coefs(:,3)=aux;

arc_middle = nrbcirc(r_P_int,[0 0],beta,pi/2);
aux=arc_middle.coefs(:,1);
arc_middle.coefs(:,1)=arc_middle.coefs(:,3);
arc_middle.coefs(:,3)=aux;
  
arc_bottom = nrbcirc(rl,[0 z_c],3*pi/2,3*pi/2+atan(a/z_c));    % Arco del Lente

arc_right = nrbcirc(r,[0 0],0,alpha);

%%----------------------------Bottom Patch-------------------------------

%Right
curve_br = line_right1;

%Top
curve_bt = arc_middle;

%Left
curve_bl = line_left1;

%Bottom
curve_bb = arc_bottom;

%Patch
Fmap_b = nrbcoons(curve_bb,curve_bt,curve_bl,curve_br);

% Calculos para graficar el Jacobiano del mapeo
[ders, ders2] = nrbderiv (Fmap_b);
npts=299;
pts = {linspace(0, 1, npts), linspace(0, 1, npts)};

% Valor de la parametrizacion evaluada en los ptos pts: matriz de orden 2 x npts^2 
[F,Jac] = geo_nurbs (Fmap_b, ders, ders2, pts, 1, 2); 

%Mean Ratio Jacobian- medida de la calidad de la parametrizacion
Zvector(:)=Jac(1,1,:).*Jac(2,2,:)-Jac(1,2,:).*Jac(2,1,:);
denom(:)=Jac(1,1,:).*Jac(1,1,:)+Jac(2,1,:).*Jac(2,1,:)+ Jac(1,2,:).*Jac(1,2,:)+Jac(2,2,:).*Jac(2,2,:); 
Jrvector=2*Zvector./denom;
ZJr_b=reshape(Jrvector,npts,npts);

norminfJac_b=max(max(max(Jac)));
MinMrJ_b =min(min(ZJr_b));

X_b=reshape(F(1,:),npts,npts);
Y_b=reshape(F(2,:),npts,npts);

%%---------------------------------Top Patch--------------------------------
%Right
curve_tr = line_right2;
%Top
curve_tt = arc_top;
%Left
curve_tl = line_left2;
%Bottom
curve_tb = arc_middle;
%Patch
Fmap_t = nrbcoons(curve_tb,curve_tt,curve_tl,curve_tr);

% Calculos para graficar el Jacobiano del mapeo
[ders, ders2] = nrbderiv (Fmap_t);
pts = {linspace(0, 1, npts), linspace(0, 1, npts)};

% Valor de la parametrizacion evaluada en los ptos pts: matriz de orden 2 x npts^2 
[F,Jac] = geo_nurbs (Fmap_t, ders, ders2, pts, 1, 2); 

%Mean Ratio Jacobian- medida de la calidad de la parametrizacion
Zvector(:)=Jac(1,1,:).*Jac(2,2,:)-Jac(1,2,:).*Jac(2,1,:);
denom(:)=Jac(1,1,:).*Jac(1,1,:)+Jac(2,1,:).*Jac(2,1,:)+ Jac(1,2,:).*Jac(1,2,:)+Jac(2,2,:).*Jac(2,2,:); 
Jrvector=2*Zvector./denom;
ZJr_t=reshape(Jrvector,npts,npts);

norminfJac_t=max(max(max(Jac)));
MinMrJ_t =min(min(ZJr_t));

X_t=reshape(F(1,:),npts,npts);
Y_t=reshape(F(2,:),npts,npts);

%%---------------------------Right Patch-------------
%Right
curve_rr = arc_right;
%Top
curve_rt = line_right2;
%Left
curve_rl = line_right1;
%Bottom
curve_rb = line_right3;
%Patch
Fmap_r = nrbcoons(curve_rb,curve_rt,curve_rl,curve_rr);

% Calculos para graficar el Jacobiano del mapeo
[ders, ders2] = nrbderiv (Fmap_r);
pts = {linspace(0, 1, npts), linspace(0, 1, npts)};

% Valor de la parametrizacion evaluada en los ptos pts: matriz de orden 2 x npts^2 
[F,Jac] = geo_nurbs (Fmap_r, ders, ders2, pts, 1, 2); 

%Mean Ratio Jacobian- medida de la calidad de la parametrizacion
Zvector(:)=Jac(1,1,:).*Jac(2,2,:)-Jac(1,2,:).*Jac(2,1,:);
denom(:)=Jac(1,1,:).*Jac(1,1,:)+Jac(2,1,:).*Jac(2,1,:)+ Jac(1,2,:).*Jac(1,2,:)+Jac(2,2,:).*Jac(2,2,:); 
Jrvector=2*Zvector./denom;
ZJr_r=reshape(Jrvector,npts,npts);

norminfJac_r=max(max(max(Jac)));
MinMrJ_r =min(min(ZJr_r));

X_r=reshape(F(1,:),npts,npts);
Y_r=reshape(F(2,:),npts,npts);

%close;
%----------------------------------------
MinMrJ=min([MinMrJ_b, MinMrJ_t, MinMrJ_r])
norminfJac = max([norminfJac_t, norminfJac_b, norminfJac_r]);
Fmap = [Fmap_b, Fmap_t, Fmap_r];

[geometry, boundaries,interfaces,~,boundary_interfaces]= mp_geo_load(Fmap);



if show_figure == 1
    figure(10), hold on
    nrbctrlplot(line_right1)
    nrbctrlplot(line_right2)
    nrbctrlplot(line_right3)
    nrbctrlplot(line_left1)
    nrbctrlplot(line_left2)
    nrbctrlplot(arc_top)
    nrbctrlplot(arc_middle)
    nrbctrlplot(arc_bottom)
    nrbctrlplot(arc_right)
    ax=gca;
    ax.FontSize=14;
    hold off

    figure(2), hold on
    nrbctrlplot(Fmap_b)
    nrbctrlplot(Fmap_t)
    nrbctrlplot(Fmap_r)
    ax=gca;
    ax.FontSize=14;
    hold off

    figure(5),hold on, grid on, view(2)
    surf (X_b, Y_b, ZJr_b,'EdgeColor','none','Facecolor','interp')
    surf (X_t, Y_t, ZJr_t,'EdgeColor','none','Facecolor','interp')
    surf (X_r, Y_r, ZJr_r,'EdgeColor','none','Facecolor','interp')    
    colorbar
    title('Mean ratio Jacobian')
    axis tight, axis equal
    ax=gca;
    ax.FontSize=14;
    hold off

end