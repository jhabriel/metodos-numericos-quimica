% EJ_06_OPTIMIZACION_SECCION_DORADA_Y_NEWTON
% Maximo de f(x) = 2 sin(x) - x^2/10 en [0, 4] (ejemplo clasico de Chapra)
% por seccion dorada (sin derivadas) y por Newton (con f' y f'').
%
% Ejecutar:  run('ejemplos/ej_06_optimizacion_seccion_dorada_y_newton.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

f   = @(x) 2*sin(x) - x.^2/10;
df  = @(x) 2*cos(x) - x/5;
d2f = @(x) -2*sin(x) - 1/5;

[xg, fg, itg] = seccion_aurea(f, 0, 4, 1e-8, 'max');
[xn, itn, hn, tipo] = newton_optimizacion(df, d2f, 2.5, 1e-10);

fprintf('Seccion dorada : x* = %.8f, f(x*) = %.8f  (%d iteraciones)\n', xg, fg, itg);
fprintf('Newton (x0=2.5): x* = %.8f, f(x*) = %.8f  (%d iteraciones) -> %s\n\n', ...
        xn, f(xn), itn, tipo);
fprintf('Valor conocido : x* = 1.42755178, f(x*) = 1.77572566\n\n');

fprintf('Newton, historia:\n   k        x_k           f''(x_k)       f''''(x_k)\n');
fprintf('%4d  %.10f  %12.3e  %10.4f\n', hn');

% Newton busca puntos estacionarios: con otro x0 puede ir a un minimo
[xm, ~, ~, tipo2] = newton_optimizacion(df, d2f, 5, 1e-10);
fprintf('\nCon x0 = 5 Newton converge a x = %.6f, que es un %s: ', xm, tipo2);
fprintf('depende del punto inicial.\n');

figure;
xx = linspace(0, 4, 300);
plot(xx, f(xx), 'b-', xg, fg, 'ro', xn, f(xn), 'kx', 'LineWidth', 1.3, 'MarkerSize', 10);
xlabel('x'); ylabel('f(x)'); grid on;
legend('f(x)', 'Seccion dorada', 'Newton', 'Location', 'south');
title('f(x) = 2 sin x - x^2/10');
