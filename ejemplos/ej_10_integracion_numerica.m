% EJ_10_INTEGRACION_NUMERICA
% Formulas cerradas de Newton-Cotes (de una sola aplicacion sobre todo [a, b]):
% trapecio, Simpson 1/3, Simpson 3/8 y Boole, para
%
%        I = integral_1^2  dV / V  = ln 2
%
% (trabajo isotermico reversible de un gas ideal: W = n R T * I)
%
% Cada formula usa n+1 nodos equiespaciados x_i = a + i h, con h = (b - a)/n:
%
%   Trapecio     (n = 1)  I ~ (h/2)   [ f0 + f1 ]
%   Simpson 1/3  (n = 2)  I ~ (h/3)   [ f0 + 4 f1 + f2 ]
%   Simpson 3/8  (n = 3)  I ~ (3h/8)  [ f0 + 3 f1 + 3 f2 + f3 ]
%   Boole        (n = 4)  I ~ (2h/45) [ 7 f0 + 32 f1 + 12 f2 + 32 f3 + 7 f4 ]
%
% Se reporta el error relativo verdadero porcentual:
%
%   e_t = | I_exacta - I_aproximada | / | I_exacta | * 100 %
%
% Ejecutar:  run('ejemplos/ej_10_integracion_numerica.m')

clear; clc; close all;

f = @(x) 1 ./ x;
a = 1; b = 2;
Iex = log(2);                      % valor exacto

%% Trapecio (n = 1)
h = (b - a) / 1;
x = a + (0:1) * h;
y = f(x);
I_trap = h / 2 * (y(1) + y(2));

%% Simpson 1/3 (n = 2)
h = (b - a) / 2;
x = a + (0:2) * h;
y = f(x);
I_s13 = h / 3 * (y(1) + 4 * y(2) + y(3));

%% Simpson 3/8 (n = 3)
h = (b - a) / 3;
x = a + (0:3) * h;
y = f(x);
I_s38 = 3 * h / 8 * (y(1) + 3 * y(2) + 3 * y(3) + y(4));

%% Boole (n = 4)
h = (b - a) / 4;
x = a + (0:4) * h;
y = f(x);
I_boole = 2 * h / 45 * (7 * y(1) + 32 * y(2) + 12 * y(3) + 32 * y(4) + 7 * y(5));

%% Resultados
I  = [I_trap, I_s13, I_s38, I_boole];
et = abs(Iex - I) / abs(Iex) * 100;          % error relativo verdadero (%)
nombres = {'Trapecio', 'Simpson 1/3', 'Simpson 3/8', 'Boole'};

fprintf('Integral exacta = ln 2 = %.10f\n\n', Iex);
fprintf('Metodo         n   Aproximacion   Error relativo verdadero (%%)\n');
for i = 1:4
    fprintf('%-12s  %2d   %.10f    %.6f\n', nombres{i}, i, I(i), et(i));
end
