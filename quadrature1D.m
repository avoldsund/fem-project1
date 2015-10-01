function I = quadrature1D(a, b, Nq, g)
% Approximates the integral of a function g(x) with integration
% limits a to b. Nq is the number of integration points.

% Get the Gaussian quadrature points x_q and the associated Gaussian
% weights w_q.
[X, W] = gauss_quadrature_1D(Nq);
I = 0;

for i = 1:Nq
    I = I + (b - a) * W(i) * g(a + (b-a) * X(i));
end

end