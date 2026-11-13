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

% Polinomio de grado 4 (5 nodos) y de grado 2 (3 nodos cercanos a 45 C)
P4 = lagrange(T, P, Tq);
P2 = lagrange(T(3:5), P(3:5), Tq);
P1 = lagrange(T(3:4), P(3:4), Tq);

fprintf('Estimacion de P(45 C):\n');
fprintf('  grado 1 (40, 50)        : %.4f kPa  (error %.4f)\n', P1, abs(P1 - Pref));
fprintf('  grado 2 (40, 50, 60)    : %.4f kPa  (error %.4f)\n', P2, abs(P2 - Pref));
fprintf('  grado 4 (todos)         : %.4f kPa  (error %.4f)\n', P4, abs(P4 - Pref));
fprintf('  valor de tablas         : %.4f kPa\n', Pref);

% --- Grafica ---
Tf = linspace(20, 60, 200);
figure; hold on;
plot(Tf, lagrange(T, P, Tf), 'b-', 'LineWidth', 1.3);
plot(T, P, 'ko', 'MarkerFaceColor', 'k');
plot(Tq, P4, 'rs', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('T (^\circC)'); ylabel('P_{vapor} (kPa)'); grid on;
legend('Polinomio de Lagrange (grado 4)', 'Datos', 'Estimacion a 45 ^\circC', 'Location', 'northwest');
title('Presion de vapor del agua');
