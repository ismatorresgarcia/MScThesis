%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%PROGRAMA PARA REPRESENTAR LA DENSIDAD DE KRIPTÓN 8+ DEL CANAL%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TODOS LOS PARÁMETROS ESTÁN EN MICRÓMETROS
r = linspace(0,60); % Variable radio del canal 
z = linspace(0,5000); % Variable longitud de propagación del haz

% PARÁMETROS DE LA PARÁBOLA DE DENSIDAD DE IONES DE KRIPTON 
r_0 = 80; % Radio máximo del canal de plasma 
r_c = 38; % Radio del canal de plasma
r_v = 90; % Radio auxiliar del canal de plasma

% REPRESENTACIÓN DE LA PARÁBOLA DE IONES DE Kr DEL CANAL DE PLASMA
N_Kr = @(r) (r<r_c).*(1+(r./r_0).^2) + ... 
    (r>=r_c).*((1+(r_c/r_0).^2).*((r_v-r)./(r_v-r_c)));
figure(2)
f = fplot(N_Kr,[0 60]);
set(f,'linewidth',1);
grid on

% LEYENDA 1
xlabel({'Radio del canal','r (\mum)'})
ylabel({'Densidad de iones de Kr','N_{Kr}'})

% PARÁMETROS DE LA PARÁBOLA DE DENSIDAD DE IONES DE KRIPTÓN 8+
z_0 = 3750; % Posición de la zona de ionización intermedia
sigmaz = 250; % Desviación estándar de la longitud de propagación
sigmar = 15; % Desviación estándar del radio del canal
[ri] = f_sigmoide(50,5,2.5,4.1); % Llamada a la rutina f_sigmoide
rLi = ri; % Guardo en memoria la sigmoide interior
% rLe = re; % Guardo en memoria la sigmoide exterior
sigmarL = 15; % Desviación estándar del ancho del canal con Kr8+
sigmaz0 = 2000; % Longitud de la primera zona de sobreionización
rL0 = 40; % Anchura máxima del canal con Kr8+
sigmar0 = 15; % Desviación estándar de r_0

% TÉRMINOS QUE APARECEN EN LA FUNCIÓN DE DENSIDAD DE Kr8+
fact_exp1i = @(r,z) exp((-0.5.*max(r,rLi(z)).^2)./sigmarL^2) ...
    ./exp(-0.5.*rLi(z).^2/sigmarL^2);
% fact_exp1e = @(r,z) exp((-0.5.*max(r,rLe(z)).^2)./sigmarL^2) ... 
%     ./exp(-0.5.*rLe(z).^2/sigmarL^2);
fact_exp2 = @(r) exp((-0.5.*max(r,rL0).^2)./sigmar0^2) ...
    ./exp(-0.5*rL0^2/sigmar0^2);
par_1 = @(r,z) 1-exp((-0.5.*(z-z_0).^2)./sigmaz^2) ... 
    .*exp(-0.5.*r.^2./sigmar^2);
par_2 = @(r,z) 1-exp((-0.5.*z.^2)./sigmaz0^2) ... 
    .*exp(-0.5.*r.^2./sigmar^2).*fact_exp2(r);

% PRODUCTOS ESCALARES INTERMEDIOS AUXILIARES
p1 = @(r,z) N_Kr(r).*par_1(r,z);
p2_i = @(r,z) fact_exp1i(r,z).*par_2(r,z);
% p2_e = @(r,z) fact_exp1e(r,z).*par_2(r,z);

% REPRESENTACIÓN 3D DE LA EVOLUCIÓN DE IONES DE Kr8+ 
N_Kr8_i = @(r,z) p1(r,z).*p2_i(r,z);
% N_Kr8_e = @(r,z) p1(r,z).*p2_e(r,z);
figure(3)
fmesh(N_Kr8_i,[0 60 0 5000])
% figure(4)
% fmesh(N_Kr8_e,[0 80 0 5000])
% zlabel({'Densidad de iones de Kr^{8+} exterior','N_{Kr^{8+}}'})

% LEYENDA
xlabel({'Radio del canal','r (\mum)'})
ylabel({'Distancia de propagación','z (\mum)'})
zlabel({'Densidad de iones de Kr^{8+} interior','N_{Kr^{8+}}'})



