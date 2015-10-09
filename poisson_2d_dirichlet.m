function u = poisson_2d_dirichlet(f, nr_of_mesh_nodes)
% The function takes a f function and a number of mesh nodes. It
% approximates the solution for u, and plots this, based on interpolation.
% The plotted solution shows u on unit disc.

% Get the triangulation elements, the vertices and edges
[p, tri, edge] = getDisk(nr_of_mesh_nodes);

% Get the stiffness matrix and the b vector before they are altered from
% the boundary conditions.
[A, b] = get_stiffness_matrix_and_b(nr_of_mesh_nodes, f, p, tri);

% Implementing Dirichlet boundary conditions
% Setting all columns and rows with an index on the edge to 0. Setting the
% diagonal elements to 1 for the zero rows. Let the elements in the
% b-vector be equal to 0 for the edge indices.
A(:, edge(:,1)) = 0;
A(edge(:,1), :) = 0;
A(edge(:,1), edge(:,1)) = eye(length(edge));
b(edge(:,1)) = 0;

% Solving the system to find the weights for u.
u = A\b;

% Plotting the solution
trimesh(tri, p(:,1), p(:,2), u)
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
trimesh(tri, p(:,1), p(:,2), error)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolated plots and solutions
% Plot of the u solution
xlin = linspace(-1, 1, nr_of_mesh_nodes);
ylin = linspace(-1, 1, nr_of_mesh_nodes);
[X, Y] = meshgrid(xlin, ylin);
f = scatteredInterpolant(p(:,1), p(:,2), u);
U = f(X,Y);
U(X.^2 + Y.^2 > 1) = 0;

figure
mesh(X,Y,U)
str=sprintf('Numerical approximation using %d points', nr_of_mesh_nodes);
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
mesh(X, Y, E)
title('Error plot')
colorbar


% Plot of u. To do this we follow a procedure for surface plots of
% nonuniformly sampled data. We create two vectors xlin and ylin that are
% uniformly spaced, and create a meshgrid based on these. We then
% interpolate the values of the function at the uniformly spaced points
% based on the function values in the old nonuniform sample points.


end