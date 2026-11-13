% EJ_05_SISTEMA_NO_LINEAL_NEWTON
% Sistema no lineal clasico de 3x3 (Burden & Faires):
%
%   3 x1 - cos(x2 x3) - 1/2                       = 0
%   x1^2 - 81 (x2 + 0.1)^2 + sin(x3) + 1.06       = 0
%   exp(-x1 x2) + 20 x3 + (10 pi - 3)/3           = 0
%
% Solucion exacta: x = [0.5; 0; -pi/6]
%
% Ejecutar:  run('ejemplos/ej_05_sistema_no_lineal_newton.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

F = @(x) [ 3*x(1) - cos(x(2)*x(3)) - 1/2;
           x(1)^2 - 81*(x(2) + 0.1)^2 + sin(x(3)) + 1.06;
           exp(-x(1)*x(2)) + 20*x(3) + (10*pi - 3)/3 ];

J = @(x) [ 3,                  x(3)*sin(x(2)*x(3)),    x(2)*sin(x(2)*x(3));
           2*x(1),             -162*(x(2) + 0.1),      cos(x(3));
           -x(2)*exp(-x(1)*x(2)), -x(1)*exp(-x(1)*x(2)), 20 ];

x0 = [0.1; 0.1; -0.1];
xex = [0.5; 0; -pi/6];

[x1, it1, h1] = newton_sistema(F, J, x0, 1e-12);   % jacobiano analitico
[x2, it2]     = newton_sistema(F, [], x0, 1e-12);  % jacobiano por dif. finitas

fprintf('Jacobiano analitico  : %d iteraciones, ||x - x_exacta|| = %.2e\n', it1, norm(x1 - xex));
fprintf('Jacobiano numerico   : %d iteraciones, ||x - x_exacta|| = %.2e\n\n', it2, norm(x2 - xex));
fprintf('x = [%s]\n', num2str(x1', ' %.10f'));
fprintf('||F(x)|| = %.2e\n\n', norm(F(x1)));

fprintf('Historia (analitico):\n   k    ||dx||       ||F||\n');
fprintf('%4d  %.3e  %.3e\n', h1');

figure;
semilogy(h1(:, 1), h1(:, 2), 'o-', h1(:, 1), h1(:, 3), 's-', 'LineWidth', 1.3);
xlabel('iteracion'); grid on; legend('||\Deltax||', '||F(x)||');
title('Newton para sistemas: convergencia cuadratica');
