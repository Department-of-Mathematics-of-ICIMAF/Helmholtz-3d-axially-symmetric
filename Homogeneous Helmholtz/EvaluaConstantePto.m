function[cxy, alphaxy, rhoxy] = EvaluaConstantePto(x, y, c_const, alpha_const, rho_const)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EvaluaConstantePto.m  -  Caso k CONSTANTE (Homogeneous Helmholtz)
% VHM octubre 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper para evaluar las magnitudes fisicas
%   c      : velocidad de propagacion del ultrasonido,
%   alpha  : coeficiente de atenuacion,
%   rho    : densidad del tejido,
% en una matriz de puntos (x, y) del dominio fisico.
%
% En el modelo con k CONSTANTE estas tres magnitudes son valores constantes
% en todo el dominio. Por ello, este helper simplemente devuelve matrices
% del mismo tamano que las entradas x e y, llenas con los valores constantes
% correspondientes.
%
% Este helper es el analogo, en el caso k constante, de la funcion
% EvaluaShoenbergPto.m del caso k variable. En el caso k variable se usa el
% spline de Shoenberg Vg(y) = sum_{i=1}^n g(tstar_i) * B_i^k(y) para
% aproximar cada magnitud constante por pedazos. Aqui, como las magnitudes
% ya son constantes, se prescinde de los splines y se devuelven valores
% constantes directamente. Mantener la misma firma de salida
% (cxy, alphaxy, rhoxy) permite que evaluaQ.m sea casi identico en ambas
% carpetas.
%
% Entrada:
%   x, y        - matrices n x m con las coordenadas fisicas de los puntos
%                 donde se evaluan las magnitudes
%   c_const     - valor constante de la velocidad de propagacion del
%                 ultrasonido (m/s)
%   alpha_const - valor constante del coeficiente de atenuacion (Np/m)
%   rho_const   - valor constante de la densidad del tejido (kg/m^3)
%
% Salida:
%   cxy         - matriz n x m con el valor c_const en todas sus entradas
%   alphaxy     - matriz n x m con el valor alpha_const en todas sus entradas
%   rhoxy       - matriz n x m con el valor rho_const en todas sus entradas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, m] = size(x);
cxy     = c_const     * ones(n, m);
alphaxy = alpha_const * ones(n, m);
rhoxy   = rho_const   * ones(n, m);
end
