% EJ_10_INTEGRACION_NUMERICA
% Formulas cerradas de Newton-Cotes (compuestas): trapecio, Simpson 1/3,
% Simpson 3/8 y Boole, para
%
%        I = integral_1^2  dV / V  = ln 2
%
% (trabajo isotermico reversible de un gas ideal: W = n R T * I)
%
% Con n subintervalos, h = (b - a)/n y nodos x_i = a + i h (i = 0, ..., n),
% cada formula es una suma ponderada de los valores y_i = f(x_i):
%
%   Trapecio     (h/2)   * sum w_i y_i,  w = [1 2 2 ... 2 1]              cualquier n
%   Simpson 1/3  (h/3)   * sum w_i y_i,  w = [1 4 2 4 2 ... 4 1]          n par
%   Simpson 3/8  (3h/8)  * sum w_i y_i,  w = [1 3 3 2 3 3 2 ... 3 3 1]    n multiplo de 3
%   Boole        (2h/45) * sum w_i y_i,  w = [7 32 12 32 14 32 12 32 14 ... 7]
%                                                                        n multiplo de 4
%
% Se usan n = 12, 24, 48, 96: multiplos de 1, 2, 3 y 4, de modo que las cuatro
% formulas son aplicables con el mismo n.
%
% Ejecutar:  run('ejemplos/ej_10_integracion_numerica.m')

clear; clc; close all;

f = @(x) 1 ./ x;
a = 1; b = 2;
Iex = log(2);

n = [12 24 48 96];
I = zeros(numel(n), 4);          % columnas: trapecio, Simpson 1/3, Simpson 3/8, Boole
for k = 1:numel(n)
    nk = n(k);
    h = (b - a) / nk;
    x = a + (0:nk) * h;          % nodos
    y = f(x);                    % valores en los nodos

    % Trapecio: pesos [1 2 2 ... 2 1]
    w = 2 * ones(1, nk + 1);
    w([1 end]) = 1;
    I(k, 1) = h / 2 * sum(w .* y);

    % Simpson 1/3: pesos [1 4 2 4 2 ... 4 1]
    w = 2 * ones(1, nk + 1);
    w(2:2:nk) = 4;
    w([1 end]) = 1;
    I(k, 2) = h / 3 * sum(w .* y);

    % Simpson 3/8: pesos [1 3 3 2 3 3 2 ... 3 3 1]
    w = 3 * ones(1, nk + 1);
    w(4:3:nk) = 2;
    w([1 end]) = 1;
    I(k, 3) = 3 * h / 8 * sum(w .* y);

    % Boole: pesos [7 32 12 32 14 32 12 32 14 ... 7]
    w = 14 * ones(1, nk + 1);
    w(2:2:nk) = 32;
    w(3:4:nk) = 12;
    w([1 end]) = 7;
    I(k, 4) = 2 * h / 45 * sum(w .* y);
end

E = abs(I - Iex);                % errores absolutos

fprintf('Integral exacta = ln 2 = %.12f\n\n', Iex);
fprintf('Valores con n = 12:\n');
fprintf('  Trapecio     : %.12f\n', I(1, 1));
fprintf('  Simpson 1/3  : %.12f\n', I(1, 2));
fprintf('  Simpson 3/8  : %.12f\n', I(1, 3));
fprintf('  Boole        : %.12f\n\n', I(1, 4));

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
