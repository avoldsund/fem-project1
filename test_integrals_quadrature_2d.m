function I = test_integrals_quadrature_2d()
% Solution of the first problem in 1 b)

fun1 = @(x, y) log(x + y);
xmin1 = 1;
xmax1 = 3;
ymin1 = @(x) 1/2 * x - 1/2;
ymax1 = @(x) x - 1;

I = integral2(fun1, xmin1, xmax1, ymin1, ymax1);

end