function [X, W] = gauss_quadrature_1D(Nq)
% Function that returns the Gaussian quadrature points and their associated
% Gaussian weights in 1D. Nq is the number of points.

X_4_point = zeros(4);
W_4_point = zeros(4);

X_4_point(1,1) = 1/2;
X_4_point(2,1) = 1/2 - sqrt(3)/6;
X_4_point(2,2) = 1/2 + sqrt(3)/6;
X_4_point(3,1) = 1/2 - sqrt(15)/10;
X_4_point(3,2) = 1/2;
X_4_point(3,3) = 1/2 + sqrt(15)/10;
X_4_point(4,1) = 1/2 - sqrt(525 + 70 * sqrt(30))/70;
X_4_point(4,2) = 1/2 - sqrt(525 - 70 * sqrt(30))/70;
X_4_point(4,3) = 1/2 + sqrt(525 - 70 * sqrt(30))/70;
X_4_point(4,4) = 1/2 + sqrt(525 + 70 * sqrt(30))/70;

W_4_point(1,1) = 1;
W_4_point(2,1) = 1/2;
W_4_point(2,2) = 1/2;
W_4_point(3,1) = 5/18;
W_4_point(3,2) = 4/9;
W_4_point(3,3) = 5/18;
W_4_point(4,1) = (18 - sqrt(30))/72;
W_4_point(4,2) = (18 + sqrt(30))/72;
W_4_point(4,3) = (18 + sqrt(30))/72;
W_4_point(4,4) = (18 - sqrt(30))/72;

X = X_4_point(Nq, 1:Nq);
W = W_4_point(Nq, 1:Nq);

end