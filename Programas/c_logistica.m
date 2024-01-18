%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%PROGRAMA PARA REPRESENTAR CURVAS LOGÍSTICAS%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = linspace(0,5); % Variable longitud de propagación del haz (mm)

% PARÁMETROS DE LA CURVA INTERIOR
ki = 50; % Tasa de decrecimiento de la curva
ri_max = 5; % Valor máximo que toma la curva
ri_min = 0; % Valor mínimo que toma la curva
z0i = 1; % El valor de z donde la curva toma su valor medio
Li = ri_max-ri_min; % Parámetro auxiliar

% PARÁMETROS DE LA CURVA EXTERIOR
ke = 2; % Tasa de decrecimiento de la curva
re_max = 3; % Valor máximo que toma la curva
re_min = 2; % Valor mínimo que toma la curva
z0e = 0; % El valor de z donde la curva toma su valor medio
Le = re_max-re_min; % Parámetro auxiliar

% FUNCIÓN RADIO DEL CANAL (SIGMOIDES)
ri = @(z) ri_min + Li./(1+exp(ki.*(z-z0i))); % Radio del canal interior  
... con iones de Kr^8+
re = @(z) re_min + Le./(1+exp(ke.*(z-z0e))); % Radio del canal exterior  
... con iones de Kr^8+

% REPRESENTACIÓN GRÁFICA DE LAS CURVAS LOGÍSTICAS
f1 = fplot(ri,[0 5]); % Frontera interior
set(f1(1),'linewidth',1);
grid on
hold on
f2 = fplot(re,[0 5]); % Frontera exterior
set(f2(1),'linewidth',1);
hold off

% LEYENDA
xlabel({'Distancia de propagación','z (mm)'})
ylabel({'Radio del canal con iones','r (\mum)'})
legend('Frontera interior','Frontera exterior')



