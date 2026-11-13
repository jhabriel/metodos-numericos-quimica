function test_ejemplos()
% Prueba de humo: cada script de ejemplos/ debe ejecutarse completo, sin
% errores y sin advertencias. Las figuras se crean ocultas y la salida de texto
% se descarta.

raiz = fileparts(fileparts(mfilename('fullpath')));
archivos = dir(fullfile(raiz, 'ejemplos', 'ej_*.m'));
assert(~isempty(archivos), 'no se encontraron ejemplos en ejemplos/');

visible = get(0, 'defaultfigurevisible');
set(0, 'defaultfigurevisible', 'off');

% Sin pantalla (por ejemplo en GitHub Actions) Octave usa el motor gnuplot y
% avisa de que esta desaconsejado. Es un aviso del entorno, no de los ejemplos.
estado_aviso = warning('off', 'Octave:gnuplot-graphics');

fallos = {};
for i = 1:numel(archivos)
    ruta = fullfile(raiz, 'ejemplos', archivos(i).name);
    lastwarn('');
    try
        ejecutar_ejemplo(ruta);
        msg = lastwarn();
        if ~isempty(msg)
            fallos{end+1} = sprintf('%s: advertencia: %s', archivos(i).name, msg); %#ok<AGROW>
        end
    catch err
        fallos{end+1} = sprintf('%s: %s', archivos(i).name, err.message); %#ok<AGROW>
    end
end

close all;
set(0, 'defaultfigurevisible', visible);
warning(estado_aviso);

if ~isempty(fallos)
    error('Fallaron %d ejemplo(s):\n  %s', numel(fallos), strjoin(fallos, sprintf('\n  ')));
end
end

function ejecutar_ejemplo(ruta)
% Los ejemplos empiezan con "clear; clc; close all;". Se ejecutan en esta
% funcion auxiliar para que ese "clear" no borre las variables de
% test_ejemplos. evalc captura (y descarta) el texto que imprimen.
evalc('run(ruta);');
end
