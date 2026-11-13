% EJ_05_SISTEMA_NO_LINEAL_NEWTON
% Sistema no lineal clasico de 3x3 (Burden & Faires):
%
%   3 x1 - cos(x2 x3) - 1/2                       = 0
%   x1^2 - 81 (x2 + 0.1)^2 + sin(x3) + 1.06       = 0
%   exp(-x1 x2) + 20 x3 + (10 pi - 3)/3           = 0
%
% Solucion exacta: x = [0.5; 0; -pi/6]
%
% Metodo de Newton para sistemas, armado en el script: en cada iteracion se
% resuelve el sistema lineal
%
%     J(x_k) dx = -F(x_k)        (con resolver_lu)
%
% y se actualiza x_{k+1} = x_k + dx, donde J es el jacobiano analitico de F.
%
% Ejecutar:  run('ejemplos/ej_05_sistema_no_lineal_newton.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

% Funcion no lineal R^3 -> R^3
F = @(x) [ 3*x(1) - cos(x(2)*x(3)) - 1/2;
           x(1)^2 - 81*(x(2) + 0.1)^2 + sin(x(3)) + 1.06;
           exp(-x(1)*x(2)) + 20*x(3) + (10*pi - 3)/3 ];

% Jacobiano R^3 -> R^(3x3):  J(i,j) = dF_i/dx_j
J = @(x) [ 3,                      x(3)*sin(x(2)*x(3)),    x(2)*sin(x(2)*x(3));
           2*x(1),                 -162*(x(2) + 0.1),      cos(x(3));
           -x(2)*exp(-x(1)*x(2)),  -x(1)*exp(-x(1)*x(2)),  20 ];

x0 = [0.1; 0.1; -0.1];  % aproximacion inicial
xex = [0.5; 0; -pi/6];  % valor exacto (de referencia)

%% Iteracion de Newton
tol = 1e-12;
maxit = 50;
x = x0;
hist = zeros(0, 3);                        % columnas: k, ||dx||, ||F(x)||
for k = 1:maxit
    dx = resolver_lu(J(x), -F(x));         % resolver J*dx = -F
    x = x + dx;
    hist(end+1, :) = [k, norm(dx), norm(F(x))];
    if norm(dx) < tol
        break
    end
end

%% Resultados
fprintf('Newton: %d iteraciones, ||x - x_exacta|| = %.2e\n\n', k, norm(x - xex));
fprintf('x = [%s]\n', num2str(x', ' %.10f'));
fprintf('||F(x)|| = %.2e\n\n', norm(F(x)));

fprintf('Historia:\n   k    ||dx||       ||F||\n');
fprintf('%4d  %.3e  %.3e\n', hist');

figure;
semilogy(hist(:, 1), hist(:, 2), 'o-', hist(:, 1), hist(:, 3), 's-', 'LineWidth', 1.3);
xlabel('iteracion'); grid on; legend('||\Deltax||', '||F(x)||');
title('Newton para sistemas: convergencia cuadratica');
