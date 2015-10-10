function I = quadLine2D(a, b, Nq, g)
% quadLine2D calculates the integral of the function g on the line segment
% going from the starting point a to the end point b. Nq is the number of
% integration points.
[X, W] = gauss_quadrature_1D(Nq);

I = 0;

for i = 1:Nq
    I = I + norm(b - a) * W(i) * g(a + (b-a) * X(i));
end

end