%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%PROGRAMA INTERFASE MATLAB-FORTRAN%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ri] = f_sigmoide(ki,ri_max,ri_min,z0i)
% function [re] = f_sigmoide(ke,re_max,re_min,z0e)

% PARÁMETROS DE LA CURVA INTERIOR
% ki: Tasa de decrecimiento de la curva interior
% ri_max: Valor máximo que toma la curva interior
% ri_min: Valor mínimo que toma la curva interior
% z0i: El valor de z donde la curva interior toma su valor medio
Li = ri_max-ri_min; % Parámetro auxiliar

% PARÁMETROS DE LA CURVA EXTERIOR
% ke: Tasa de decrecimiento de la curva
% re_max: Valor máximo que toma la curva exterior
% re_min: Valor mínimo que toma la curva interior
% z0e: El valor de z donde la curva exterior toma su valor medio
% Le = re_max-re_min; % Parámetro auxiliar

% FUNCIÓN RADIO DEL CANAL (SIGMOIDES)
ri = @(z) ri_min + Li./(1+exp(ki.*(z-z0i))); % Radio del canal interior
... con iones de Kr^8+
% re = @(z) re_min + Le./(1+exp(ke.*(z-z0e))); % Radio del canal exterior 
... con iones de Kr^8+

% REPRESENTACIÓN GRÁFICA DE LAS CURVAS LOGÍSTICAS
figure(1)
fig1 = fplot(ri,[0 5]); % Frontera interior 
set(fig1(1),'linewidth',1);
grid on
% fig2 = fplot(re,[0 5]); % Frontera exterior
% set(fig2(1),'linewidth',1);
% hold off

% LEYENDA
xlabel({'Distancia de propagación','z (mm)'})
ylabel({'Radio del canal','r (\mum)'})
legend('Frontera interior')

end

