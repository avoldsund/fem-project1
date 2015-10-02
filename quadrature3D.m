function I = quadrature3D(p1, p2, p3, p4, Nq, g)
% Approximates the integral of a function g(x1, x2, x3) of the domain Omega, made
% up of the tetrahedron p1, p2, p3 and p4. Nq is the number of Gaussian integration
% points.

% Get the 3D Gauss quadrature weights and integration points.
[Zeta, W] = gauss_quadrature_3D(Nq);
B = [p1(1) - p4(1), p2(1) - p4(1), p3(1) - p4(1);
    p1(2) - p4(2), p2(2) - p4(2), p3(2) - p4(2);
    p1(3) - p4(3), p2(3) - p4(3), p3(3) - p4(3)];
Jacobian = abs(det(B));
I = 0;

for q = 1:Nq
    % Do the change of variable from zeta_1, zeta_2 and zeta_3 to x.
    x = Zeta(q, 1) * p1 + Zeta(q, 2) * p2 + Zeta(q, 3) * p3 + Zeta(q, 4) * p4;
    I = I + W(q) * g(x) * Jacobian;
end

end