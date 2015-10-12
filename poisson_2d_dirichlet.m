function poisson_2d_dirichlet(nr_of_mesh_nodes)
% The function takes a number of mesh nodes. It
% approximates the solution for u, and plots this, based on interpolation.
% The plotted solution shows u on unit disc.

f = @(x) -8 * pi * cos(2 * pi * (x(1)^2 + x(2)^2)) + 16 * pi^2 * (x(1)^2 + x(2)^2) * sin(2 * pi * (x(1)^2 + x(2)^2));

% Get the triangulation elements, the vertices and edges
[p, tri, edge] = getDisk(nr_of_mesh_nodes);

% Get the stiffness matrix and the b vector before they are altered from
% the boundary conditions.
[A, b] = get_stiffness_matrix_and_load_vector_2D(nr_of_mesh_nodes, f, p, tri);

% Implementing Dirichlet boundary conditions
% Setting all columns and rows with an index on the edge to 0. Setting the
% diagonal elements to 1 for the zero rows. Let the elements in the
% b-vector be equal to 0 for the edge indices.
boundary = edge(:,1);
A(boundary, :) = 0;
A(boundary, boundary) = speye(length(edge));
b(boundary) = 0;

% Solving the system to find the weights for u.
u = A\b;

% Plotting the solution
trimesh(tri, p(:,1), p(:,2), full(u))
str = sprintf('Numerical approximation using %d points', nr_of_mesh_nodes);
title(str)

u_analytical = zeros(nr_of_mesh_nodes, 1);
k = @(x) sin(2 * pi * (x(1)^2 + x(2)^2));
for i = 1 : nr_of_mesh_nodes
    u_analytical(i) = k(p(i,:));
end
    
error = abs(u - u_analytical);
max_error = max(error);
% Plot of the error
figure
trimesh(tri, p(:,1), p(:,2), full(error))
colorbar
str = sprintf('Error plot. Max error: %f ', max_error);
title(str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolated plots and solutions
% Plot of the u solution
xlin = linspace(-1, 1, nr_of_mesh_nodes);
ylin = linspace(-1, 1, nr_of_mesh_nodes);
[X, Y] = meshgrid(xlin, ylin);
f = scatteredInterpolant(p(:,1), p(:,2), full(u));
U = f(X,Y);
U(X.^2 + Y.^2 > 1) = 0;

figure
mesh(X,Y,U)
str=sprintf('Numerical approximation using %d points and interpolation', nr_of_mesh_nodes);
title(str)

% Plot of analytical solution
analytical_solution = @(X, Y) sin(2 * pi * (X.^2 + Y.^2));
U_analytical = analytical_solution(X, Y);
U_analytical(X.^2 + Y.^2 > 1) = 0;

figure
mesh(X, Y, U_analytical)
str=sprintf('Analytical solution');
title(str)

% Plot of absolute error in the mesh
f = scatteredInterpolant(p(:,1), p(:,2), error);
E = f(X,Y);
E(X.^2 + Y.^2 > 1) = 0;
figure
mesh(X, Y, E)
title('Error plot for interpolated solution and error')
colorbar

end