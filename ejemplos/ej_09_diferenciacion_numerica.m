% EJ_09_DIFERENCIACION_NUMERICA
% Primera y segunda derivada de f(x) = exp(x) sin(x) en x0 = 1 con diferencias
% finitas, y estudio del error al reducir el paso h.
%
%   f'(x)  = exp(x) (sin x + cos x)
%   f''(x) = 2 exp(x) cos x
%
% Ejecutar:  run('ejemplos/ej_09_diferenciacion_numerica.m')

clear; clc; close all;
addpath(genpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src')));

f = @(x) exp(x) .* sin(x);
x0 = 1;
d1 = exp(x0) * (sin(x0) + cos(x0));
d2 = 2 * exp(x0) * cos(x0);

h = 0.4 ./ 2.^(0:6);
ead = zeros(size(h)); eat = ead; ece = ead; e2 = ead;
for k = 1:numel(h)
    ead(k) = abs(dif_primera(f, x0, h(k), 'adelante') - d1);
    eat(k) = abs(dif_primera(f, x0, h(k), 'atras')    - d1);
    ece(k) = abs(dif_primera(f, x0, h(k), 'central')  - d1);
    e2(k)  = abs(dif_segunda_central(f, x0, h(k))     - d2);
end

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
