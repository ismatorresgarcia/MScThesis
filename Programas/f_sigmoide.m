%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%PROGRAMA INTERFASE MATLAB-FORTRAN%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ri] = f_sigmoide(ki,ri_max,ri_min,z0i) % Parámetros sigmoide interior
function [re] = f_sigmoide(ke,re_max,re_min,z0e) % Parámetros sigmoide exterior

% PARÁMETROS DE LA CURVA INTERIOR Y EXTERIOR
Li = ri_max-ri_min; Le = re_max-re_min; % Parámetros auxiliares

% FUNCIÓN RADIO DEL CANAL (SIGMOIDES)
ri = @(z) ri_min + Li./(1+exp(ki.*(z-z0i))); % Radio del canal interior
... con iones de Kr^8+
re = @(z) re_min + Le./(1+exp(ke.*(z-z0e))); % Radio del canal exterior 
... con iones de Kr^8+

% REPRESENTACIÓN GRÁFICA DE LAS CURVAS LOGÍSTICAS
figure(1)
fig1 = fplot(ri,[0 5]); % Frontera interior 
set(fig1(1),'linewidth',1);
grid on
fig2 = fplot(re,[0 5]); % Frontera exterior
set(fig2(1),'linewidth',1);

% LEYENDA
xlabel({'Distancia de propagación','z (mm)'})
ylabel({'Radio del canal','r (\mum)'})
legend('Frontera interior')

