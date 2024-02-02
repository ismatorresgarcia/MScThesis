%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%PROGRAMA PARA REPRESENTAR LA EXPONENCIAL A TROZOS%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = linspace(0,60); % Variable radio del canal (en micrómetros)

% PARÁMETROS DE LA FUNCIÓN EXPONENCIAL
r_L = 15;     % Radio del canal con Kr8+ para r<r_u (en micrómetros)
r_u = 18;     % Radio del canal con Kr8+ para r>=r_u (en micrómetros)
sig_L = 15;   % Desviación estándar para r_L (en micrómetros)
sig_u = 17;   % Desviación estándar para r_u (en micrómetros)

% FUNCIÓN A TROZOS DEL RADIO DEL CANAL
f_L = @(r) (exp(-0.5*((max(r,r_L).^2)/sig_L^2)))/ ...
    (exp(-0.5*(r_L/sig_L)^2));                       % Trozo para r<r_u
A = f_L(r_u)*exp(0.5*(r_u/sig_u)^2);                 % Condición de continuidad
f_u = @(r) A*exp(-0.5*(r/sig_u).^2);                 % Trozo para r>=r_u
f = @(r) (r<r_u).*f_L(r) + (r>=r_u).*f_u(r);         % Función exponencial

% REPRESENTACIÓN GRÁFICA DE LA FUNCIÓN A TROZOS
fig = fplot(f); % Gráfica de la función a trozos
xlim([0 60])
ylim([0 1.1])
set(fig(1),'linewidth',1,'color','r');
grid on
% LEYENDA
xlabel({'Radio del canal','r (\mum)'})
ylabel({'f(r)','Adimensional'})



