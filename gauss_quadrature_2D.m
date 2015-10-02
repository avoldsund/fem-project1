function [Zeta, W] = gauss_quadrature_2D(Nq)
% Function that returns the Gaussian quadrature points and their associated
% Gaussian weights in 2D. Nq is the number of points.

switch Nq
    case 1
        W = 1/2;
        Zeta = [1/3 1/3 1/3];
    case 3
        W = [1/6 1/6 1/6];
        Zeta = [1/2 1/2 0; 1/2 0 1/2; 0 1/2 1/2];
    case 4
        W = [-9/32 25/96 25/96 25/96];
        Zeta = [1/3 1/3 1/3; 3/5 1/5 1/5; 1/5 3/5 1/5; 1/5 1/5 3/5];
    otherwise
        error('Number of integration points can only equal 1, 3 or 4.')
        
end