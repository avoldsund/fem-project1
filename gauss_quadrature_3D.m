function [Zeta, W] = gauss_quadrature_3D(Nq)
% Function that returns the Gaussian quadrature points and their associated
% Gaussian weights in 3D. Nq is the number of points.
a = 1/4 + 3 * sqrt(5)/20;
b = 1/4 - sqrt(5)/20;

switch Nq
    case 1
        W = 1/6;
        Zeta = [1/4 1/4 1/4 1/4];
    case 4
        W = [1/24 1/24 1/24 1/24];
        Zeta = [a b b b; b a b b; b b a b; b b b a];
    case 5
        W = [-4/30 9/120 9/120 9/120 9/120];
        Zeta = [1/4 1/4 1/4 1/4; 1/2 1/6 1/6 1/6; 1/6 1/2 1/6 1/6; 1/6 1/6 1/2 1/6; 1/6 1/6 1/6 1/2];
    otherwise
        error('Number of integration points can only equal 1, 4 or 5.')
        
end