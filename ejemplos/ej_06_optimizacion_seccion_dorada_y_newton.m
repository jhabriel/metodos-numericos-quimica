% EJ_06_OPTIMIZACION_SECCION_DORADA_Y_NEWTON
% Maximo de f(x) = 2 sin(x) - x^2/10 en [0, 4] (ejemplo clasico de Chapra)
% por seccion dorada (sin derivadas) y por Newton (con g' y g'').
%
% Los metodos de optimizacion de este repositorio MINIMIZAN. Para hallar un
% maximo de f se minimiza g(x) = -f(x); el valor maximo de f es -g(x*).
%
% Newton para optimizar es el mismo metodo de Newton para raices, aplicado a
% g'(x) = 0:   x_{k+1} = x_k - g'(x_k) / g''(x_k)   ->  newton(dg, d2g, x0).
% El signo de g'' en el punto hallado dice si es minimo (g'' > 0) o maximo (g'' < 0).
%
% Ejecutar:  run('ejemplos/ej_06_optimizacion_seccion_dorada_y_newton.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

f = @(x) 2*sin(x) - x.^2/10;               % funcion a MAXIMIZAR

% Se minimiza g = -f
g   = @(x) -f(x);
dg  = @(x) -2*cos(x) + x/5;
d2g = @(x) 2*sin(x) + 1/5;

[xg, gg, itg] = seccion_aurea(g, 0, 4, 1e-8);
[xn, itn, hn] = newton(dg, d2g, 2.5, 1e-10);     % raiz de g'

fprintf('Seccion dorada : x* = %.8f, f(x*) = %.8f  (%d iteraciones)\n', xg, -gg, itg);
fprintf('Newton (x0=2.5): x* = %.8f, f(x*) = %.8f  (%d iteraciones)\n', xn, f(xn), itn);
fprintf('   g''''(x*) = %.4f > 0: minimo de g, es decir, maximo de f\n\n', d2g(xn));
fprintf('Valor conocido : x* = 1.42755178, f(x*) = 1.77572566\n\n');

fprintf('Newton, historia:\n   k        x_k           g''(x_k)\n');
fprintf('%4d  %.10f  %12.3e\n', hn');

% Newton busca puntos estacionarios: con otro x0 puede ir a un maximo de g
xm = newton(dg, d2g, 5, 1e-10);
fprintf('\nCon x0 = 5 Newton converge a x = %.6f, donde g''''(x) = %.4f', xm, d2g(xm));
if d2g(xm) > 0
    fprintf(' > 0 (minimo de g).\n');
else
    fprintf(' < 0 (maximo de g, NO es lo buscado):\n');
    fprintf('hay que revisar el signo de g'''' y el punto inicial.\n');
end

figure;
xx = linspace(0, 4, 300);
plot(xx, f(xx), 'b-', xg, -gg, 'ro', xn, f(xn), 'kx', 'LineWidth', 1.3, 'MarkerSize', 10);
xlabel('x'); ylabel('f(x)'); grid on;
legend('f(x)', 'Seccion dorada', 'Newton', 'Location', 'south');
title('f(x) = 2 sin x - x^2/10');
