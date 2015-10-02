function I = quadrature2D(p1, p2, p3, Nq, g)
% Approximates the integral of a function g(x1, x2) of the domain Omega, made
% up of the points p1, p2 and p3. Nq is the number of Gaussian integration
% points.

% Get the 2D Gauss quadrature weights and integration points.
[Zeta, W] = gauss_quadrature_2D(Nq);
Jacobian = abs((p1(1) - p3(1)) * (p2(2) - p3(2)) - (p2(1) - p3(1)) * (p1(2) - p3(2)));
I = 0;

for q = 1:Nq
    % Do the change of variable from zeta_1, zeta_2 and zeta_3 to x.
    x = Zeta(q, 1) * p1 + Zeta(q, 2) * p2 + Zeta(q, 3) * p3;
    I = I + W(q) * g(x) * Jacobian;
end

end