function test_integracion()
% Trapecio, Simpson 1/3, Simpson 3/8 y Boole.

f = @(x) sin(x); a = 0; b = pi; ex = 2;

assert(abs(trapecio(f, a, b, 100) - ex) < 2e-4, 'trapecio');
assert(abs(simpson13(f, a, b, 20) - ex) < 1e-5, 'simpson13');
assert(abs(simpson38(f, a, b, 21) - ex) < 1e-4, 'simpson38');
assert(abs(boole(f, a, b, 16) - ex) < 1e-5, 'boole');

% Exactitud: S1/3 y S3/8 exactas hasta grado 3; Boole hasta grado 5
c3 = @(x) x.^3 + 2 * x.^2 - x + 1;  I3 = 1/4 + 2/3 - 1/2 + 1;   % en [0,1]
assert(abs(simpson13(c3, 0, 1, 2) - I3) < 1e-14, 'simpson13 cubica');
assert(abs(simpson38(c3, 0, 1, 3) - I3) < 1e-14, 'simpson38 cubica');
c5 = @(x) x.^5; assert(abs(boole(c5, 0, 1, 4) - 1/6) < 1e-14, 'boole grado 5');

% Restricciones sobre n
for k = 1:3
    ok = false;
    try
        switch k
            case 1, simpson13(f, a, b, 5);
            case 2, simpson38(f, a, b, 5);
            case 3, boole(f, a, b, 6);
        end
    catch
        ok = true;
    end
    assert(ok, 'n invalido debe fallar');
end
end
