% EJ_08_AJUSTE_DE_CURVAS
% Regresion lineal por minimos cuadrados, armando las ECUACIONES NORMALES.
%
% Para un modelo lineal en los parametros,
%
%     y ~ a1*phi1(x) + a2*phi2(x) + ... + am*phim(x),
%
% se construye la matriz Z (una columna por funcion base, una fila por dato)
% y se minimiza ||y - Z a||^2, lo que lleva al sistema lineal
%
%     (Z' * Z) a = Z' * y
%
% que se resuelve con resolver_lu. Se muestran tres casos:
%   A) Regresion lineal simple   : curva de calibracion (ley de Beer-Lambert)
%   B) Regresion polinomica      : Cp(T) de un gas (valores aproximados CO2)
%   C) Regresion lineal multiple : rendimiento vs temperatura y concentracion
%
% Requisito: haber ejecutado `iniciar` al abrir MATLAB/Octave (ver README).
% Ejecutar:  run('ejemplos/ej_08_ajuste_de_curvas.m')

clear; clc; close all;

% Coeficiente de determinacion: R^2 = 1 - Sr/St
%   Sr = suma de cuadrados de los residuos, St = suma de cuadrados respecto a la media
R2 = @(y, yaj) 1 - sum((y - yaj).^2) / sum((y - mean(y)).^2);

%% A) Recta: Absorbancia = a1 + a2 * concentracion
c   = [0 2 4 6 8 10]';                         % mg/L
Abs = [0.003 0.152 0.301 0.449 0.604 0.748]';

Z = [ones(size(c)), c];                        % columnas: 1, c
a = resolver_lu(Z' * Z, Z' * Abs);             % ecuaciones normales
r2A = R2(Abs, Z * a);

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

%% B) Polinomios: Cp(T)
T  = [300 400 500 600 700 800 900 1000]';      % K
Cp = [37.2 41.3 44.6 47.3 49.6 51.4 52.9 54.3]';   % J/(mol K)

Z2 = [ones(size(T)), T, T.^2];                 % columnas: 1, T, T^2
Z3 = [ones(size(T)), T, T.^2, T.^3];           % columnas: 1, T, T^2, T^3
a2 = resolver_lu(Z2' * Z2, Z2' * Cp);
a3 = resolver_lu(Z3' * Z3, Z3' * Cp);
r2_2 = R2(Cp, Z2 * a2);
r2_3 = R2(Cp, Z3 * a3);

fprintf('B) Cp(T) (coeficientes en potencias ascendentes de T):\n');
fprintf('   grado 2: coef = [%s],  R^2 = %.6f\n', num2str(a2', ' %.4e'), r2_2);
fprintf('   grado 3: coef = [%s],  R^2 = %.6f\n\n', num2str(a3', ' %.4e'), r2_3);

figure;
Tg = linspace(300, 1000, 200)';
plot(T, Cp, 'ko', 'MarkerFaceColor', 'k'); hold on;
plot(Tg, a2(1) + a2(2)*Tg + a2(3)*Tg.^2, 'b-', 'LineWidth', 1.3);
plot(Tg, a3(1) + a3(2)*Tg + a3(3)*Tg.^2 + a3(4)*Tg.^3, 'r--', 'LineWidth', 1.3);
xlabel('T (K)'); ylabel('C_p (J/mol K)'); grid on;
legend('Datos', sprintf('Grado 2 (R^2=%.5f)', r2_2), sprintf('Grado 3 (R^2=%.5f)', r2_3), ...
       'Location', 'southeast');
title('B) Regresion polinomica');

%% C) Multiple: rendimiento = a1 + a2*T + a3*C
Temp = [60 60 60 70 70 70 80 80 80 90 90 90]';        % C
Conc = [ 1  2  3  1  2  3  1  2  3  1  2  3]';        % mol/L
Rend = [38.4 41.4 46.3 45.8 50.5 53.6 54.1 58.6 61.5 61.7 66.2 70.4]';   % %

Z = [ones(size(Temp)), Temp, Conc];            % columnas: 1, T, C
aM = resolver_lu(Z' * Z, Z' * Rend);
yM = Z * aM;
r2M = R2(Rend, yM);

fprintf('C) Rend = %.3f + %.4f*T + %.4f*C,   R^2 = %.5f\n', aM(1), aM(2), aM(3), r2M);
Rpred = aM(1) + aM(2) * 75 + aM(3) * 2.5;
fprintf('   Prediccion para T = 75 C, C = 2.5 mol/L: %.2f %%\n', Rpred);

figure;

% Plano ajustado en el espacio (T, C, Rend), con los residuos
subplot(1, 2, 1); hold on;
[TT, CC] = meshgrid(linspace(55, 95, 15), linspace(0.5, 3.5, 15));
ZZ = aM(1) + aM(2) * TT + aM(3) * CC;
mesh(TT, CC, ZZ, 'EdgeColor', [0.4 0.6 0.9]);
hidden off;                 % malla transparente: se ven puntos y residuos detras
for i = 1:numel(Rend)
    plot3([Temp(i) Temp(i)], [Conc(i) Conc(i)], [Rend(i) yM(i)], 'r-', 'LineWidth', 1.5);   % residuo
end
plot3(Temp, Conc, Rend, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);
xlabel('T (^\circC)'); ylabel('C (mol/L)'); zlabel('Rendimiento (%)');
grid on; view(-35, 25);
title('C) Plano ajustado y residuos');

% Observado vs ajustado
subplot(1, 2, 2); hold on;
plot(Rend, yM, 'ko', 'MarkerFaceColor', 'k');
lim = [min(Rend) max(Rend)];
plot(lim, lim, 'r-', 'LineWidth', 1.2);
xlabel('Rendimiento observado (%)'); ylabel('Rendimiento ajustado (%)'); grid on;
axis equal; title(sprintf('Observado vs ajustado (R^2 = %.4f)', r2M));
