% EJ_07_INTERPOLACION_CON_LAGRANGE
% Presion de vapor del agua (kPa) a distintas temperaturas (C). Se interpola
% con el polinomio de Lagrange para estimar la presion a 45 C.
% Valor de tablas a 45 C: ~9.59 kPa.
%
% Ejecutar:  run('ejemplos/ej_07_interpolacion_con_lagrange.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

T = [20 30 40 50 60];                      % C
P = [2.339 4.246 7.384 12.352 19.946];     % kPa
Tq = 45;
Pref = 9.593;

% Coeficientes del polinomio de Lagrange (potencias descendentes, como polyval)
% Se centra y escala la temperatura (t = (T - 40)/10) para evitar mal
% condicionamiento al expandir el polinomio.
tf = @(T) (T - 40) / 10;
[p4, L4] = lagrange(tf(T), P);                 % grado 4 (5 nodos)
p2 = lagrange(tf(T(3:5)), P(3:5));             % grado 2 (40, 50, 60)
p1 = lagrange(tf(T(3:4)), P(3:4));             % grado 1 (40, 50)

P4 = polyval(p4, tf(Tq));
P2 = polyval(p2, tf(Tq));
P1 = polyval(p1, tf(Tq));

fprintf('Coeficientes (grado 4, en t = (T-40)/10, potencias descendentes):\n');
fprintf('  [%s]\n\n', num2str(p4, ' %.5f'));

fprintf('Estimacion de P(45 C):\n');
fprintf('  grado 1 (40, 50)        : %.4f kPa  (error %.4f)\n', P1, abs(P1 - Pref));
fprintf('  grado 2 (40, 50, 60)    : %.4f kPa  (error %.4f)\n', P2, abs(P2 - Pref));
fprintf('  grado 4 (todos)         : %.4f kPa  (error %.4f)\n', P4, abs(P4 - Pref));
fprintf('  valor de tablas         : %.4f kPa\n', Pref);

% --- Graficas ---
Tg = linspace(20, 60, 200);
figure;
subplot(1, 2, 1); hold on;
plot(Tg, polyval(p4, tf(Tg)), 'b-', 'LineWidth', 1.3);
plot(T, P, 'ko', 'MarkerFaceColor', 'k');
plot(Tq, P4, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('T (^\circC)'); ylabel('P_{vapor} (kPa)'); grid on;
legend('Polinomio de Lagrange (grado 4)', 'Datos', 'Estimacion a 45 ^\circC', 'Location', 'northwest');
title('Presion de vapor del agua');

subplot(1, 2, 2); hold on;
for i = 1:numel(T)
    plot(Tg, polyval(L4(i, :), tf(Tg)), 'LineWidth', 1.2);
end
plot(T, 0 * T, 'k.', 'MarkerSize', 12);
xlabel('T (^\circC)'); ylabel('L_i(T)'); grid on;
title('Polinomios base de Lagrange');
