% EJ_10_INTEGRACION_NUMERICA
% Formulas cerradas de Newton-Cotes (compuestas): trapecio, Simpson 1/3,
% Simpson 3/8 y Boole, para
%
%        I = integral_1^2  dV / V  = ln 2
%
% (trabajo isotermico reversible de un gas ideal: W = n R T * I)
%
% Se usan n = 12, 24, 48, 96 subintervalos: multiplos de 1, 2, 3 y 4, de modo
% que las cuatro formulas son aplicables con el mismo n.
%
% Ejecutar:  run('ejemplos/ej_10_integracion_numerica.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

f = @(x) 1 ./ x;
a = 1; b = 2;
Iex = log(2);

n = [12 24 48 96];
E = zeros(numel(n), 4);
for k = 1:numel(n)
    E(k, 1) = abs(trapecio(f, a, b, n(k)) - Iex);
    E(k, 2) = abs(simpson13(f, a, b, n(k)) - Iex);
    E(k, 3) = abs(simpson38(f, a, b, n(k)) - Iex);
    E(k, 4) = abs(boole(f, a, b, n(k)) - Iex);
end

fprintf('Integral exacta = ln 2 = %.12f\n\n', Iex);
fprintf('Valores con n = 12:\n');
fprintf('  Trapecio     : %.12f\n', trapecio(f, a, b, 12));
fprintf('  Simpson 1/3  : %.12f\n', simpson13(f, a, b, 12));
fprintf('  Simpson 3/8  : %.12f\n', simpson38(f, a, b, 12));
fprintf('  Boole        : %.12f\n\n', boole(f, a, b, 12));

fprintf('Errores absolutos:\n');
fprintf('   n     trapecio    Simpson1/3  Simpson3/8     Boole\n');
fprintf('%4d  %11.3e %11.3e %11.3e %11.3e\n', [n', E]');

fprintf('\nOrden observado (n -> 2n): ');
ord = mean(log2(E(1:end-1, :) ./ E(2:end, :)));
fprintf('%.2f  ', ord); fprintf('\n(teorico: 2, 4, 4, 6)\n');

h = (b - a) ./ n;
figure;
loglog(h, E(:, 1), 'o-', h, E(:, 2), 's-', h, E(:, 3), '^-', h, E(:, 4), 'd-', 'LineWidth', 1.3);
xlabel('h'); ylabel('error absoluto'); grid on;
legend('Trapecio', 'Simpson 1/3', 'Simpson 3/8', 'Boole', 'Location', 'southeast');
title('Newton-Cotes cerrado: error vs h');
