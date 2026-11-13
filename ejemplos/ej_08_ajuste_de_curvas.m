% EJ_08_AJUSTE_DE_CURVAS
% Regresion lineal generalizada (minimos cuadrados) en tres casos:
%   A) Regresion lineal simple   : curva de calibracion (ley de Beer-Lambert)
%   B) Regresion polinomica      : Cp(T) de un gas (valores aproximados CO2)
%   C) Regresion lineal multiple : rendimiento vs temperatura y concentracion
%
% Ejecutar:  run('ejemplos/ej_08_ajuste_de_curvas.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

uno = @(x) ones(size(x, 1), 1);

%% A) Recta: Absorbancia = a0 + a1 * concentracion
c   = [0 2 4 6 8 10]';                         % mg/L
Abs = [0.003 0.152 0.301 0.449 0.604 0.748]';
[a, r2A, yA] = regresion_lineal_general(c, Abs, {uno, @(x) x});

fprintf('A) Recta: A = %.5f + %.5f c,   R^2 = %.5f\n', a(1), a(2), r2A);
Amuestra = 0.380;
fprintf('   Muestra con A = %.3f  ->  c = %.3f mg/L\n\n', Amuestra, (Amuestra - a(1)) / a(2));

figure;
cc = linspace(0, 10, 100)';
plot(c, Abs, 'ko', 'MarkerFaceColor', 'k'); hold on;
plot(cc, a(1) + a(2) * cc, 'b-', 'LineWidth', 1.3);
xlabel('c (mg/L)'); ylabel('Absorbancia'); grid on;
legend('Datos', sprintf('Ajuste (R^2 = %.4f)', r2A), 'Location', 'northwest');
title('A) Curva de calibracion');

%% B) Polinomios: Cp(T)   (T escalada a kK para un sistema bien condicionado)
T  = [300 400 500 600 700 800 900 1000]';      % K
Cp = [37.2 41.3 44.6 47.3 49.6 51.4 52.9 54.3]';   % J/(mol K)
t  = T / 1000;

base2 = {uno, @(x) x, @(x) x.^2};
base3 = {uno, @(x) x, @(x) x.^2, @(x) x.^3};
[a2, r2_2] = regresion_lineal_general(t, Cp, base2);
[a3, r2_3] = regresion_lineal_general(t, Cp, base3);

fprintf('B) Cp(T) con t = T/1000:\n');
fprintf('   grado 2: coef = [%s],  R^2 = %.6f\n', num2str(a2', ' %.3f'), r2_2);
fprintf('   grado 3: coef = [%s],  R^2 = %.6f\n\n', num2str(a3', ' %.3f'), r2_3);

figure;
tt = linspace(0.3, 1.0, 200)';
plot(T, Cp, 'ko', 'MarkerFaceColor', 'k'); hold on;
plot(1000 * tt, a2(1) + a2(2)*tt + a2(3)*tt.^2, 'b-', 'LineWidth', 1.3);
plot(1000 * tt, a3(1) + a3(2)*tt + a3(3)*tt.^2 + a3(4)*tt.^3, 'r--', 'LineWidth', 1.3);
xlabel('T (K)'); ylabel('C_p (J/mol K)'); grid on;
legend('Datos', sprintf('Grado 2 (R^2=%.5f)', r2_2), sprintf('Grado 3 (R^2=%.5f)', r2_3), ...
       'Location', 'southeast');
title('B) Regresion polinomica');

%% C) Multiple: rendimiento = a0 + a1*T + a2*C
Temp = [60 60 60 70 70 70 80 80 80 90 90 90]';        % C
Conc = [ 1  2  3  1  2  3  1  2  3  1  2  3]';        % mol/L
Rend = [38.4 41.4 46.3 45.8 50.5 53.6 54.1 58.6 61.5 61.7 66.2 70.4]';   % %

X = [Temp, Conc];
baseM = {uno, @(X) X(:, 1), @(X) X(:, 2)};
[aM, r2M, yM] = regresion_lineal_general(X, Rend, baseM);

fprintf('C) Rend = %.3f + %.4f*T + %.4f*C,   R^2 = %.5f\n', aM(1), aM(2), aM(3), r2M);
Rpred = aM(1) + aM(2) * 75 + aM(3) * 2.5;
fprintf('   Prediccion para T = 75 C, C = 2.5 mol/L: %.2f %%\n', Rpred);

figure;
plot(Rend, yM, 'ko', 'MarkerFaceColor', 'k'); hold on;
lim = [min(Rend) max(Rend)];
plot(lim, lim, 'r-', 'LineWidth', 1.2);
xlabel('Rendimiento observado (%)'); ylabel('Rendimiento ajustado (%)'); grid on;
axis equal; title(sprintf('C) Regresion multiple (R^2 = %.4f)', r2M));
