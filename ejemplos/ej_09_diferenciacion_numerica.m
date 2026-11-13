% EJ_09_DIFERENCIACION_NUMERICA
% Primera y segunda derivada de f(x) = exp(x) sin(x) en x0 = 1 con diferencias
% finitas, y estudio del error al reducir el paso h.
%
%   f'(x)  = exp(x) (sin x + cos x)
%   f''(x) = 2 exp(x) cos x
%
% Formulas (h = paso):
%   adelante   f'(x) ~ ( f(x+h) - f(x)   ) / h             error O(h)
%   atras      f'(x) ~ ( f(x)   - f(x-h) ) / h             error O(h)
%   central    f'(x) ~ ( f(x+h) - f(x-h) ) / (2h)          error O(h^2)
%   2da central f''(x) ~ ( f(x+h) - 2 f(x) + f(x-h) ) / h^2  error O(h^2)
%
% Ejecutar:  run('ejemplos/ej_09_diferenciacion_numerica.m')

clear; clc; close all;

f = @(x) exp(x) .* sin(x);
x0 = 1;
d1 = exp(x0) * (sin(x0) + cos(x0));        % derivada exacta
d2 = 2 * exp(x0) * cos(x0);                % segunda derivada exacta

h = 0.4 ./ 2.^(0:6);                       % pasos: 0.4, 0.2, 0.1, ...

% Las formulas se evaluan para todos los h a la vez (h es un vector)
dadl = (f(x0 + h) - f(x0)) ./ h;                       % adelante
datr = (f(x0) - f(x0 - h)) ./ h;                       % atras
dcen = (f(x0 + h) - f(x0 - h)) ./ (2 * h);             % central
d2c  = (f(x0 + h) - 2 * f(x0) + f(x0 - h)) ./ h.^2;    % segunda derivada, central

% Errores absolutos
ead = abs(dadl - d1);
eat = abs(datr - d1);
ece = abs(dcen - d1);
e2  = abs(d2c - d2);

fprintf('f''(1) exacta = %.8f,   f''''(1) exacta = %.8f\n\n', d1, d2);
fprintf('Errores absolutos:\n');
fprintf('     h        adelante     atras       central     2da central\n');
fprintf('%8.5f  %11.3e %11.3e %11.3e %11.3e\n', [h; ead; eat; ece; e2]);

% Orden observado: log2(e(h) / e(h/2))
orden = @(e) mean(log2(e(1:end-1) ./ e(2:end)));
fprintf('\nOrden observado:  adelante %.2f,  atras %.2f,  central %.2f,  2da central %.2f\n', ...
        orden(ead), orden(eat), orden(ece), orden(e2));
fprintf('(teorico: 1, 1, 2, 2)\n');

figure;
loglog(h, ead, 'o-', h, eat, 's-', h, ece, '^-', h, e2, 'd-', 'LineWidth', 1.3);
xlabel('h'); ylabel('error absoluto'); grid on;
legend('f'' adelante', 'f'' atras', 'f'' central', 'f'''' central', 'Location', 'southeast');
title('Error de la diferenciacion numerica');
