function test_iniciar()
% iniciar.m debe dejar las funciones de src/ disponibles en el path.

raiz = fileparts(fileparts(mfilename('fullpath')));

% Se parte de una sesion "limpia": sin src/ en el path
rmpath(genpath(fullfile(raiz, 'src')));
assert(exist('biseccion', 'file') == 0, 'precondicion: biseccion no debe estar en el path');

evalc('run(fullfile(raiz, ''iniciar.m''));');      % evalc descarta el mensaje de iniciar

% Funciones de distintas carpetas de src/
funciones = {'biseccion', 'resolver_lu', 'seccion_aurea', 'lagrange', 'euler_implicito'};
for i = 1:numel(funciones)
    assert(exist(funciones{i}, 'file') == 2, ['iniciar no dejo disponible ' funciones{i}]);
end
end
