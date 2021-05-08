function B = odot( b)
    % SE2ODOT( b) returns B such that B = b\odot\xi = \xi\wedge b for SE2
    % elements.
    B = [so2alg.wedge(1) * b(1:2), b(3)*eye(2);zeros(1,3)];
end