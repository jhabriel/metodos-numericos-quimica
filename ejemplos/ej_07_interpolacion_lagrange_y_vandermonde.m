% EJ_07_INTERPOLACION_LAGRANGE_Y_VANDERMONDE
% Presion de vapor del agua (kPa) a distintas temperaturas (C). Se construye el
% polinomio interpolante de dos maneras y se comparan:
%
%   1) Matriz de Vandermonde: se imponen P(T_i) = P_i y se resuelve V a = P.
%   2) Polinomio de Lagrange: P(T) = sum_i P_i L_i(T), expandido a coeficientes.
%
% El polinomio interpolante de n puntos es unico, asi que ambos metodos deben
% dar los mismos coeficientes (salvo diferencias de redondeo). Se estima la
% presion a 45 C (valor de tablas: ~9.59 kPa).
%
% Ejecutar:  run('ejemplos/ej_07_interpolacion_lagrange_y_vandermonde.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

T = [20 30 40 50 60]';                     % C
P = [2.339 4.246 7.384 12.352 19.946]';    % kPa
Tq = 45;
Pref = 9.593;
n = numel(T);

%% 1) Vandermonde:  P(T) = a1 T^4 + a2 T^3 + a3 T^2 + a4 T + a5
% Condiciones P(T_i) = P_i  ->  sistema lineal V a = P, una fila por nodo
V = zeros(n, n);
for j = 1:n
    V(:, j) = T .^ (n - j);                % columnas: T^4, T^3, T^2, T, 1
end
aV = resolver_lu(V, P)';                   % fila, igual que lagrange

disp('Matriz de Vandermonde V (columnas: T^4, T^3, T^2, T, 1):');
disp(V)

%% 2) Lagrange (coeficientes y polinomios base)
[aL, L] = lagrange(T, P);

%% Comparacion de los coeficientes
fprintf('Coeficientes (potencias descendentes):\n');
fprintf('            Vandermonde         Lagrange        dif. relativa\n');
for j = 1:n
    fprintf('  a%d  %18.10e  %18.10e   %9.2e\n', j, aV(j), aL(j), abs(aV(j) - aL(j)) / abs(aL(j)));
end
fprintf('\n||a_Vandermonde - a_Lagrange|| / ||a_Lagrange|| = %.2e\n', norm(aV - aL) / norm(aL));
fprintf('(el polinomio es unico: la diferencia es solo error de redondeo)\n\n');

%% Comparacion: error de cada polinomio en los nodos (deberia ser 0)
fprintf('max |P(T_i) - P_i| en los nodos:\n');
fprintf('  Vandermonde : %.2e\n', max(abs(polyval(aV, T) - P)));
fprintf('  Lagrange    : %.2e\n\n', max(abs(polyval(aL, T) - P)));

%% Estimaciones a 45 C con distintos grados (nodos mas cercanos a 45 C)
fprintf('Estimacion de P(45 C), valor de tablas = %.4f kPa:\n', Pref);
fprintf('  grado   nodos            Vandermonde   Lagrange     error\n');
subconjuntos = {3:4, 3:5, 1:5};            % grado 1, 2 y 4
for s = 1:numel(subconjuntos)
    idx = subconjuntos{s};
    m = numel(idx);
    Tm = T(idx);
    Vm = zeros(m, m);
    for j = 1:m
        Vm(:, j) = Tm .^ (m - j);
    end
    pV = resolver_lu(Vm, P(idx))';
    pL = lagrange(Tm, P(idx));
    PV = polyval(pV, Tq);
    PL = polyval(pL, Tq);
    fprintf('  %d       %-14s   %.6f     %.6f     %.4f\n', m - 1, ...
            sprintf('%d-%d C', T(idx(1)), T(idx(end))), PV, PL, abs(PV - Pref));
end

%% Graficas
Tg = linspace(20, 60, 200);
figure;
subplot(1, 2, 1); hold on;
plot(Tg, polyval(aV, Tg), 'b-', 'LineWidth', 1.5);
plot(Tg, polyval(aL, Tg), 'r--', 'LineWidth', 1.2);
plot(T, P, 'ko', 'MarkerFaceColor', 'k');
plot(Tq, polyval(aV, Tq), 'ms', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('T (^\circC)'); ylabel('P_{vapor} (kPa)'); grid on;
legend('Vandermonde', 'Lagrange', 'Datos', 'Estimacion a 45 ^\circC', 'Location', 'northwest');
title('Presion de vapor del agua (las curvas coinciden)');

subplot(1, 2, 2); hold on;
for i = 1:n
    plot(Tg, polyval(L(i, :), Tg), 'LineWidth', 1.2);
end
plot(T, 0 * T, 'k.', 'MarkerSize', 12);
xlabel('T (^\circC)'); ylabel('L_i(T)'); grid on;
title('Polinomios base de Lagrange');
