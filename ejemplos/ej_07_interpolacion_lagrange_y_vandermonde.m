% EJ_07_INTERPOLACION_LAGRANGE_Y_VANDERMONDE
% Presion de vapor del agua (kPa) a distintas temperaturas (C). Se construye el
% polinomio interpolante de dos maneras y se estima la presion a 45 C
% (valor de tablas: ~9.59 kPa):
%
%   1) Matriz de Vandermonde: se imponen P(T_i) = P_i y se resuelve V a = P.
%   2) Polinomio de Lagrange: con la funcion lagrange.
%
% En ambos casos se obtienen los coeficientes del polinomio en potencias
% descendentes, que se evaluan con polyval.
%
% Requisito: haber ejecutado `iniciar` al abrir MATLAB/Octave (ver README).
% Ejecutar:  run('ejemplos/ej_07_interpolacion_lagrange_y_vandermonde.m')

clear; clc; close all;

T = [20 30 40 50 60]';                     % C
P = [2.339 4.246 7.384 12.352 19.946]';    % kPa
Tq = 45;
Pref = 9.593;

%% 1) Vandermonde:  P(T) = a1 T^4 + a2 T^3 + a3 T^2 + a4 T + a5
% Condiciones P(T_i) = P_i  ->  sistema lineal V a = P, una fila por nodo
V = [T.^4, T.^3, T.^2, T, ones(size(T))];
aV = resolver_lu(V, P)';                   % fila, igual que lagrange

disp('Matriz de Vandermonde V (columnas: T^4, T^3, T^2, T, 1):');
disp(V)
fprintf('Coeficientes (Vandermonde): [%s]\n\n', num2str(aV, ' %.6e'));

%% 2) Lagrange: coeficientes y polinomios base
[aL, L] = lagrange(T, P);
fprintf('Coeficientes (Lagrange)   : [%s]\n\n', num2str(aL, ' %.6e'));

%% Estimacion de la presion a 45 C
PV = polyval(aV, Tq);
PL = polyval(aL, Tq);
fprintf('P(45 C) con Vandermonde : %.4f kPa  (error %.4f)\n', PV, abs(PV - Pref));
fprintf('P(45 C) con Lagrange    : %.4f kPa  (error %.4f)\n', PL, abs(PL - Pref));
fprintf('Valor de tablas         : %.4f kPa\n', Pref);

%% Graficas
Tg = linspace(20, 60, 200);
figure;
subplot(1, 2, 1); hold on;
plot(Tg, polyval(aV, Tg), 'b-', 'LineWidth', 1.5);
plot(Tg, polyval(aL, Tg), 'r--', 'LineWidth', 1.2);
plot(T, P, 'ko', 'MarkerFaceColor', 'k');
plot(Tq, PV, 'ms', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('T (^\circC)'); ylabel('P_{vapor} (kPa)'); grid on;
legend('Vandermonde', 'Lagrange', 'Datos', 'Estimacion a 45 ^\circC', 'Location', 'northwest');
title('Presion de vapor del agua');

% Polinomios base de Lagrange: L_i(T_i) = 1 y L_i(T_j) = 0 para j ~= i
subplot(1, 2, 2); hold on;
for i = 1:numel(T)
    plot(Tg, polyval(L(i, :), Tg), 'LineWidth', 1.2);
end
plot(T, 0 * T, 'k.', 'MarkerSize', 12);
xlabel('T (^\circC)'); ylabel('L_i(T)'); grid on;
title('Polinomios base de Lagrange');
