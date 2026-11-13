% EJ_11_EDO_SIMPLE
% Reaccion de primer orden:   dC/dt = -k C,   C(0) = C0
% Solucion analitica:         C(t) = C0 exp(-k t)
% Se resuelve con Euler explicito, Euler implicito y Heun.
% (Euler implicito se usa para sistemas lineales y' = A y; aqui A = -k.)
%
% Ejecutar:  run('ejemplos/ej_11_edo_simple.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

k = 0.5;             % 1/s
C0 = 2;              % mol/L
tf = 10;             % s
f = @(t, C) -k * C;  % funcion anonima de la var. indep. y de la var. dep.
A = -k;              % dC/dt = A*C (lo usa Euler implicito)
Cex = @(t) C0 * exp(-k * t);

%% Una corrida con n = 20 pasos
n = 20;
[t, Ce] = euler_explicito(f, [0 tf], C0, n);
[~, Ci] = euler_implicito(A, [0 tf], C0, n);
[~, Ch] = heun(f, [0 tf], C0, n);

fprintf('n = %d pasos (h = %.2f s)\n', n, tf / n);
fprintf('  t     exacta   E.explicito  E.implicito   Heun\n');
idx = 1:4:n+1;
fprintf('%5.1f  %8.5f  %9.5f  %9.5f  %9.5f\n', [t(idx), Cex(t(idx)), Ce(idx), Ci(idx), Ch(idx)]');

%% Estudio de convergencia: error en t = tf
ns = [10 20 40 80 160];
err = zeros(numel(ns), 3);
for j = 1:numel(ns)
    [~, A1] = euler_explicito(f, [0 tf], C0, ns(j));
    [~, A2] = euler_implicito(A, [0 tf], C0, ns(j));
    [~, A3] = heun(f, [0 tf], C0, ns(j));
    err(j, :) = abs([A1(end), A2(end), A3(end)] - Cex(tf));
end
fprintf('\nError en t = %g s:\n   n    E.explicito  E.implicito     Heun\n', tf);
fprintf('%4d  %11.3e %11.3e %11.3e\n', [ns', err]');
fprintf('Orden observado (ultimo refinamiento): %.2f  %.2f  %.2f   (teorico 1, 1, 2)\n', ...
        log2(err(end-1, :) ./ err(end, :)));

%% Graficas
figure;
subplot(1, 2, 1);
tt = linspace(0, tf, 200);
plot(tt, Cex(tt), 'k-', t, Ce, 'ro--', t, Ci, 'gs--', t, Ch, 'b^--', 'LineWidth', 1.2);
xlabel('t (s)'); ylabel('C (mol/L)'); grid on;
legend('Exacta', 'E. explicito', 'E. implicito', 'Heun');
title('A \rightarrow productos (1er orden), n = 20');

subplot(1, 2, 2);
loglog(tf ./ ns, err(:, 1), 'o-', tf ./ ns, err(:, 2), 's-', tf ./ ns, err(:, 3), '^-', 'LineWidth', 1.2);
xlabel('h'); ylabel('error en t = t_f'); grid on;
legend('E. explicito', 'E. implicito', 'Heun', 'Location', 'southeast');
title('Convergencia');
