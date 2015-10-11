function I = quadLine2D(a, b, Nq, g)
% quadLine2D calculates the integral of the function g on the line segment
% going from the starting point a to the end point b. Nq is the number of
% integration points.

x = @(s) a + s * (b - a);
h = @(s) g(x(s));
I = norm(b - a) * quadrature1D(0, 1, Nq, h);

end