% Grafico de la solucion de Helmholtz 3D con simetria axial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparacion IgA vs FEM para tejido x capas
figure(100), hold on, grid on
load('absprHelmholtz3DLoidel_n1e100m1e300n3e400_tejidoxcapasNdof189201.mat','Z','abspr');
plot(Z,abspr/5e5,'-b','LineWidth',2),

% Z=[];
% abspr=[];
% load('absprHelmholtz3DVHM_n1e100m1e300n3e400_tejidoxcapasNdof189201.mat','Z','abspr')
% plot(Z,abspr,'-r','LineWidth',2),

COMSOLpr=importdata('COMSOLHelmholtz3D_abs_presionlambdaentre7.txt');
plot(COMSOLpr(:,1),COMSOLpr(:,2),'-k','LineWidth',2)

xlabel('Z','FontName','TimesNewRoman','FontSize',14,'FontWeight','bold')
ylabel('abs(p)','FontName','TimesNewRoman','FontSize',14,'FontWeight','bold')
ax=gca; ax.FontSize=12;
title('Abs. value of acoustic presure 3D model multilayer tissue ')
legend('IgA Loidel Ndof=189 201', 'FEM Moreno Ndof=1 600 000')

%legend('IgA Loidel Ndof=189 201','IgA Vicky Ndof=189 201', 'FEM Moreno Ndof=1 600 000')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparacion IgA vs FEM para tejido phantom
figure(200), hold on, grid on

load('absprHelmholtz3DVHM_n1e100m1e300n3e400_tejidophantomNdof189201.mat','Z','abspr')
plot(Z,abspr,'-r','LineWidth',2),

load('absprHelmholtz3DVHM_n1e120m1e360n3e480_tejidophantomNdof272641.mat','Z','abspr')
plot(Z,abspr,'-b','LineWidth',2),


COMSOLpr_phantom=importdata('COMSOLphantom8c58.txt');
plot(COMSOLpr_phantom(:,1),COMSOLpr_phantom(:,2),'-k','LineWidth',2)

xlabel('Z','FontName','TimesNewRoman','FontSize',14,'FontWeight','bold')
ylabel('abs(p)','FontName','TimesNewRoman','FontSize',14,'FontWeight','bold')
ax=gca; ax.FontSize=12;
title('Abs. value of acoustic presure 3D model phantom tissue ')
legend('IgA Vicky Ndof=189 201','IgA Vicky Ndof=272 641', 'FEM Moreno Ndof=1 618 539')
%legend('IgA Vicky Ndof=189 201','FEM Moreno Ndof=1 618 539')




