clear all
close all
%% Inicializar variables y parámetros

%Parámetros
Pe  =   5;
Pe1 =   5;      %Numero de peclet para transferencia de masa
Pe2 =   5;      %Numero de peclet para transferencia de calor
beta=   1.5;    %Coeficiente adimensional de transferencia de calor
b   =   12;     %Aumento adiabatico adimensional de la temperatura
Gama=   20;     %Energía de activación adimensional
Da  =   0.09;    %Damko¨hler number, aproximación, en el articulo se calcula

% Discretización espacial y temporal
z = linspace(0, 1, 100);   % Dominio en z (0 a 1) con 100 puntos
dz = z(2) - z(1);
dt = 0.01*dz;  % Paso de tiempo 
t_max = 1;  % Tiempo máximo
num_steps = round(t_max / dt);  % Número de pasos de tiempo ajustado
t = linspace(0, t_max, num_steps);   % Dominio temporal ajustado

% Condiciones iniciales para x1 y x2 (pueden ajustarse)
x1 = zeros(length(z), length(t));   % Convección
x2 = zeros(length(z), length(t));   % Temperatura adimensional
x1(:, 1) = 0;  % Condición inicial para x1
x2(:, 1) = 0;  % Condición inicial para x2

% Matrices para almacenar derivadas
dx1dt = zeros(length(z), length(t));
dx2dt = zeros(length(z), length(t));
dx1dz = zeros(length(z), length(t));
dx2dz = zeros(length(z), length(t));
d2x1dz2 = zeros(length(z), length(t));
d2x2dz2 = zeros(length(z), length(t));

% Matrices para derivadas fraccionarias
% alpha
alpha = [0.5 0.8 0.9 1 0.3];

% Nuestra
%Comprobar si es que está correcto sugerencia
x1_own = zeros(length(z), length(t), length(alpha));   % Concentración
x2_own = zeros(length(z), length(t), length(alpha));   % Temperatura

dfx1dz_own = zeros(length(z), length(t), length(alpha));
dfx2dz_own = zeros(length(z), length(t), length(alpha));
df2x1dz2_own = zeros(length(z), length(t), length(alpha)); 
df2x2dz2_own = zeros(length(z), length(t), length(alpha));
dfx2dt_own = zeros(length(z), length(t), length(alpha));
dfx1dt_own = zeros(length(z), length(t), length(alpha));

% Khalil
x1_kha = zeros(length(z), length(t), length(alpha));   % Concentración
x2_kha = zeros(length(z), length(t), length(alpha));   % Temperatura

dfx1dz_kha = zeros(length(z), length(t), length(alpha));
dfx2dz_kha = zeros(length(z), length(t), length(alpha));
df2x1dz2_kha = zeros(length(z), length(t), length(alpha)); 
df2x2dz2_kha = zeros(length(z), length(t), length(alpha));
dfx2dt_kha = zeros(length(z), length(t), length(alpha));
dfx1dt_kha = zeros(length(z), length(t), length(alpha));
%% Generar la función inicial
for n = 1:length(t)-1
        for i = 2:length(z)-1
            % Derivadas espaciales por diferencias finitas
            dx1dz(i,n) = (x1(i+1, n) - x1(i-1, n)) / (2*dz);%Central diference
            dx2dz(i,n) = (x2(i+1, n) - x2(i-1, n)) / (2*dz);
            d2x1dz2(i,n) = (x1(i+1, n) - 2*x1(i, n) + x1(i-1, n)) / dz^2;
            d2x2dz2(i,n) = (x2(i+1, n) - 2*x2(i, n) + x2(i-1, n)) / dz^2;
            
                    % Implementación de Runge-Kutta 4 para x1 y x2
            % Paso 1: Calcular k1
            k1_x1 = dt * (Pe1^-1 * d2x1dz2(i,n) - dx1dz(i,n) + Da * (1 - x1(i, n)) * exp(x2(i, n) / (1 + x2(i, n) / Gama)));
            k1_x2 = dt * (Pe2^-1 * d2x2dz2(i,n) - dx2dz(i,n) - beta * x2(i, n) + b * Da * (1 - x1(i, n)) * exp(x2(i, n) / (1 + x2(i, n) / Gama)));
    
            % Paso 2: Calcular k2
            k2_x1 = dt * (Pe1^-1 * d2x1dz2(i,n) - dx1dz(i,n) + Da * (1 - (x1(i, n) + k1_x1/2)) * exp((x2(i, n) + k1_x2/2) / (1 + (x2(i, n) + k1_x2/2) / Gama)));
            k2_x2 = dt * (Pe2^-1 * d2x2dz2(i,n) - dx2dz(i,n) - beta * (x2(i, n) + k1_x2/2) + b * Da * (1 - (x1(i, n) + k1_x1/2)) * exp((x2(i, n) + k1_x2/2) / (1 + (x2(i, n) + k1_x2/2) / Gama)));
    
            % Paso 3: Calcular k3
            k3_x1 = dt * (Pe1^-1 * d2x1dz2(i,n) - dx1dz(i,n) + Da * (1 - (x1(i, n) + k2_x1/2)) * exp((x2(i, n) + k2_x2/2) / (1 + (x2(i, n) + k2_x2/2) / Gama)));
            k3_x2 = dt * (Pe2^-1 * d2x2dz2(i,n) - dx2dz(i,n) - beta * (x2(i, n) + k2_x2/2) + b * Da * (1 - (x1(i, n) + k2_x1/2)) * exp((x2(i, n) + k2_x2/2) / (1 + (x2(i, n) + k2_x2/2) / Gama)));
    
            % Paso 4: Calcular k4
            k4_x1 = dt * (Pe1^-1 * d2x1dz2(i,n) - dx1dz(i,n) + Da * (1 - (x1(i, n) + k3_x1)) * exp((x2(i, n) + k3_x2) / (1 + (x2(i, n) + k3_x2) / Gama)));
            k4_x2 = dt * (Pe2^-1 * d2x2dz2(i,n) - dx2dz(i,n) - beta * (x2(i, n) + k3_x2) + b * Da * (1 - (x1(i, n) + k3_x1)) * exp((x2(i, n) + k3_x2) / (1 + (x2(i, n) + k3_x2) / Gama)));
    
            % Actualización de x1 y x2 usando el promedio ponderado de los incrementos
            x1(i, n+1) = x1(i, n) + (k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1) / 6;
            x2(i, n+1) = x2(i, n) + (k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2) / 6;
    
            dx1dt(i,n) = (x1(i, n+1) - x1(i, n)) / (dt);
            dx2dt(i,n) = (x2(i, n+1) - x2(i, n)) / (dt);
            % Verificar si las variables están tomando valores inestables
            if abs(x1(i, n+1)) > 1e5 || abs(x2(i, n+1)) > 1e5
                fprintf('Inestabilidad detectada en el paso de tiempo %d\n', n);
                return;
            end 
        end
   % Condiciones de borde
   dx1dz(1) = Pe1 * x1(1, n); dx2dz(1) = Pe2 * x2(1, n);  % En z=0
   dx1dz(end) = 0; dx2dz(end) = 0;  % En z=1
end


%% 
fx1_RL=ones(9702, 1);
fx2_RL=ones(9702, 1);
x1_RL = zeros(98, 99, length(alpha));
x2_RL = zeros(98, 99, length(alpha));

h = waitbar(0, 'Procesando...'); % Crear la barra de progreso
totalSteps = length(alpha) * 100; % Total de iteraciones
stepCounter = 0; % Contador de iteraciones realizadas

for a = 1:length(alpha)
    for i = 1:100 
        % Actualizar contador de pasos e incrementar la barra de progreso
        stepCounter = stepCounter + 1;
        waitbar(stepCounter / totalSteps, h, sprintf('Procesando... %.2f%%', (stepCounter / totalSteps) * 100));
        x1_RL(:,:,a) = solveFractionalSystemCRL(alpha(a), alpha(a) + 1, 100, 100, 0.01, dz, fx1_RL(:));
        x2_RL(:,:,a) = solveFractionalSystemTRL(alpha(a), alpha(a) + 1, 100, 100, 0.01, dz, fx2_RL(:));
        fx1_RL = Da.*(1 + x1_RL(:,:,a)) .* exp(x2_RL(:,:,a) ./ (1 + x2_RL(:,:,a) ./ Gama));
        fx2_RL = -beta .* x2_RL(:,:,a) + b * Da .* (1 + x1_RL(:,:,a)) .* exp(x2_RL(:,:,a) ./ (1 + x2_RL(:,:,a) ./ Gama));
    end
end

close(h); % Cerrar la barra de progreso
%% alpha 0.3
fx1_RL03=ones(9702, 1);
fx2_RL03=ones(9702, 1);
x1_RL03 = zeros(98, 99);
x2_RL03 = zeros(98, 99);

for i = 1:100 
        % Actualizar contador de pasos e incrementar la barra de progreso
        x1_RL03(:,:) = solveFractionalSystemCRL(0.3, 1.3, 100, 100, 0.01, dz, fx1_RL03(:));
        x2_RL03(:,:) = solveFractionalSystemTRL(0.3, 1.3, 100, 100, 0.01, dz, fx2_RL03(:));
        fx1_RL03 = Da.*(1 + x1_RL03(:,:)) .* exp(x2_RL03(:,:) ./ (1 + x2_RL03(:,:) ./ Gama));
        fx2_RL03 = -beta .* x2_RL03(:,:) + b * Da .* (1 + x1_RL03(:,:)) .* exp(x2_RL03(:,:) ./ (1 + x2_RL03(:,:) ./ Gama));
end

%%
fx1=ones(9702, 1);
fx2=ones(9702, 1);
x1_M = zeros(98, 99);
x2_M = zeros(98, 99);

for i = 1:100 
        % Actualizar contador de pasos e incrementar la barra de progreso
        x1_M(:,:) = solveFractionalSystemCRL(1, 2, 100, 100, 0.01, dz, fx1(:));
        x2_M(:,:) = solveFractionalSystemTRL(1, 2, 100, 100, 0.01, dz, fx2(:));
        fx1 = Da.*(1 + x1_M(:,:)) .* exp(x2_M(:,:) ./ (1 + x2_M(:,:) ./ Gama));
        fx2 = -beta .* x2_M(:,:) + b * Da .* (1 + x1_M(:,:)) .* exp(x2_M(:,:) ./ (1 + x2_M(:,:) ./ Gama));
end


%% Interpolación para graficos de podbluni
x1_pf = zeros(length(z), length(t), length(alpha));
x2_pf = zeros(length(z), length(t), length(alpha));
[Z_Rl,T_Rl] = meshgrid(linspace(0,1,99),linspace(0,1,98));
[Z,T] = meshgrid(linspace(0,1,num_steps),linspace(0,1,100));
[Zo,To] = meshgrid(linspace(0,1,100),linspace(0,1,98));
for a = 1:length(alpha)-1
    x1_pf(:,:,a) = interp2(Z_Rl, T_Rl, x1_RL(:,:,a), Z, T, 'linear');
    x2_pf(:,:,a) = interp2(Z_Rl, T_Rl, x2_RL(:,:,a), Z, T, 'linear');
end
x1_M = interp2(Z_Rl, T_Rl, x1_M(:,:), Z, T, 'linear');
x2_M = interp2(Z_Rl, T_Rl, x2_M(:,:), Z, T, 'linear');
%%
x1_pf03 = interp2(Z_Rl, T_Rl, x1_RL03(:,:), Z, T, 'linear');
x2_pf03 = interp2(Z_Rl, T_Rl, x2_RL03(:,:), Z, T, 'linear');

%% Resolver la FPDR de nuestra definición
tic
for n = 1:length(t)-1
    for i= 2:length(z)-1
        for a = 1:length(alpha)
            %Nuestra
            dfx1dz_own(i,n,a) = (x1_own(i+1, n,a) - x1_own(i-1, n,a)) / (2*dz);%Central diference
            dfx2dz_own(i,n,a) = (x2_own(i+1, n,a) - x2_own(i-1, n,a)) / (2*dz);
            df2x1dz2_own(i,n,a) = (x1_own(i+1, n,a) - 2*x1_own(i, n,a) + x1_own(i-1, n,a)) / dz^2;
            df2x2dz2_own(i,n,a) = (x2_own(i+1, n,a) - 2*x2_own(i, n,a) + x2_own(i-1, n,a)) / dz^2;
    
            x1_own(i, n+1,a) = x1_own(i, n,a) + dt/exp(t(n)*(1-alpha(a))) * (Pe1^-1 *exp(z(i)*(1-alpha(a)))* df2x1dz2_own(i,n,a) - exp(z(i)*(1-alpha(a)))*dfx1dz_own(i,n,a) + Da * (1 - x1_own(i, n,a)) * exp(x2_own(i, n,a) / (1 + x2_own(i, n,a) / Gama)));
            x2_own(i, n+1,a) = x2_own(i, n,a) + dt/exp(t(n)*(1-alpha(a))) * (Pe2^-1 *exp(z(i)*(1-alpha(a)))* df2x2dz2_own(i,n,a) - exp(z(i)*(1-alpha(a)))*dfx2dz_own(i,n,a) - beta * x2_own(i, n,a) + b * Da * (1 - x1_own(i, n,a)) * exp(x2_own(i, n,a) / (1 + x2_own(i, n,a) / Gama)));
     
            dfx1dt_own(i,n,a) = (x1_own(i, n+1,a) - x1_own(i, n,a)) / (dt);
            dfx2dt_own(i,n,a) = (x2_own(i, n+1,a) - x2_own(i, n,a)) / (dt);
        end
    end
end
toc
n_finalt = length(t);
n_finalz = length(z);
%% khalil, no funciona
for n = 1:length(t)-1
    for i= 2:length(z)-1
        for a = 1:length(alpha)
            %Nuestra
            dfx1dz_kha(i,n,a) = (x1_kha(i+1, n,a) - x1_kha(i-1, n,a)) / (2*dz);%Central diference
            dfx2dz_kha(i,n,a) = (x2_kha(i+1, n,a) - x2_kha(i-1, n,a)) / (2*dz);
            df2x1dz2_kha(i,n,a) = (x1_kha(i+1, n,a) - 2*x1_kha(i, n,a) + x1_kha(i-1, n,a)) / dz^2;
            df2x2dz2_kha(i,n,a) = (x2_kha(i+1, n,a) - 2*x2_kha(i, n,a) + x2_kha(i-1, n,a)) / dz^2;
    
            x1_kha(i, n+1,a) = x1_kha(i, n,a) + dt/(max(t(i), 1e-6)^(1 - alpha(a))) * (Pe1^-1 *(max(z(i), 1e-6)^(1 - alpha(a)))* df2x1dz2_kha(i,n,a) - (max(z(i), 1e-6)^(1 - alpha(a)))*dfx1dz_kha(i,n,a) + Da * (1 - x1_kha(i, n,a)) * exp(x2_kha(i, n,a) / (1 + x2_kha(i, n,a) / Gama)));
            x2_kha(i, n+1,a) = x2_kha(i, n,a) + dt/(max(t(i), 1e-6)^(1 - alpha(a))) * (Pe2^-1 *(max(z(i), 1e-6)^(1 - alpha(a)))* df2x2dz2_kha(i,n,a) - (max(z(i), 1e-6)^(1 - alpha(a)))*dfx2dz_kha(i,n,a) - beta * x2_kha(i, n,a) + b * Da * (1 - x1_kha(i, n,a)) * exp(x2_kha(i, n,a) / (1 + x2_kha(i, n,a) / Gama)));
     
            dfx1dt_kha(i,n,a) = (x1_kha(i, n+1,a) - x1_kha(i, n,a)) / (dt);
            dfx2dt_kha(i,n,a) = (x2_kha(i, n+1,a) - x2_kha(i, n,a)) / (dt);
            if isnan(x1_kha(i, n+1,a)) || isnan(x2_kha(i, n+1,a))
                fprintf('NaN detectado en alfa=%.2f, n=%d, i=%d\n', alpha(a), n, i);
                return;
            end

        end
    end
end
%% Guardar datos para Texstudio
podbluni_05 = [z',x1_pf(:, n_finalt-1,1)];
podbluni_08 = [z',x1_pf(:, n_finalt-1,2)];
podbluni_09 = [z',x1_pf(:, n_finalt-1,3)];
nuestra_05 = [z',x1_own(:, n_finalt-1,1)];
nuestra_08 = [z',x1_own(:, n_finalt-1,2)];
nuestra_09 = [z',x1_own(:, n_finalt-1,3)];
nuestra_1 = [z',x1_own(:, n_finalt-1,4)];
podbluni_1 = [z',x1_M(:,n_finalt-1)];
save('data_podbluni_alfa05.dat', 'podbluni_05', '-ascii', '-double', '-tabs');
save('data_podbluni_alfa08.dat', 'podbluni_08', '-ascii', '-double', '-tabs');
save('data_podbluni_alfa09.dat', 'podbluni_09', '-ascii', '-double', '-tabs');
save('data_podbluni_alfa1.dat', 'podbluni_1', '-ascii', '-double', '-tabs');
save('data_nuestra_alfa05.dat', 'nuestra_05', '-ascii', '-double', '-tabs');
save('data_nuestra_alfa08.dat', 'nuestra_08', '-ascii', '-double', '-tabs');
save('data_nuestra_alfa09.dat', 'nuestra_09', '-ascii', '-double', '-tabs');
save('data_nuestra_alfa1.dat', 'nuestra_1', '-ascii', '-double', '-tabs');
%%
x1_euler = [z', x1(:, n_finalt -1)];
save('dataX1.dat', 'x1_euler', '-ascii', '-double', '-tabs');
%%
podbluni2_05 = [z',x2_pf(:, n_finalt-1,1)];
podbluni2_08 = [z',x2_pf(:, n_finalt-1,2)];
podbluni2_09 = [z',x2_pf(:, n_finalt-1,3)];
nuestra2_05 = [z',x2_own(:, n_finalt-1,1)];
nuestra2_08 = [z',x2_own(:, n_finalt-1,2)];
nuestra2_09 = [z',x2_own(:, n_finalt-1,3)];
nuestra2_1 = [z',x2_own(:, n_finalt-1,4)];
podbluni2_1 = [z',x2_M(:,n_finalt-1)];
x2_euler = [z', x2(:, n_finalt -1)];
save('dataX2.dat', 'x2_euler', '-ascii', '-double', '-tabs');
save('data_podbluni_alfa205.dat', 'podbluni2_05', '-ascii', '-double', '-tabs');
save('data_podbluni_alfa208.dat', 'podbluni2_08', '-ascii', '-double', '-tabs');
save('data_podbluni_alf2a09.dat', 'podbluni2_09', '-ascii', '-double', '-tabs');
save('data_podbluni_alfa21.dat', 'podbluni2_1', '-ascii', '-double', '-tabs');
save('data_nuestra_alfa205.dat', 'nuestra2_05', '-ascii', '-double', '-tabs');
save('data_nuestra_alfa208.dat', 'nuestra2_08', '-ascii', '-double', '-tabs');
save('data_nuestra_alfa209.dat', 'nuestra2_09', '-ascii', '-double', '-tabs');
save('data_nuestra_alfa21.dat', 'nuestra2_1', '-ascii', '-double', '-tabs');
%%
podbluni_03 = [z',x1_pf03(:, n_finalt-1)];
nuestra_03 = [z',x1_own(:, n_finalt-1,5)];
save('datapodblunialfa03.dat', 'podbluni_03', '-ascii', '-double', '-tabs');
save('datanuestraalfa03.dat', 'nuestra_03', '-ascii', '-double', '-tabs');
%% Gráficos
% alpha = [0.5 0.8 0.9 1];
figure;
hold on 
plot(z, x1_pf(:, n_finalt-1,1), '-g', 'LineWidth', 1.5, 'DisplayName', 'Podlubny \alpha = 0.9');
% plot(z, x1_M(:, n_finalt-1), '--b', 'LineWidth', 1.5, 'DisplayName', 'Podlubny');
plot(z, x1_own(:, n_finalt-1,1), '-c', 'LineWidth', 1.5, 'DisplayName', 'Euler \alpha = 0.9');
% plot(z, x1(:, n_finalt-1), '--k', 'LineWidth', 1.5, 'DisplayName', 'Euler');
hold off
xlabel('Z');
ylabel('C');
title('Conversion dimentionless form at t = 1');
legend('Podlubny \alpha = 0.5', 'Euler \alpha = 0.5', 'Euler khalil \alpha = 0.5');
