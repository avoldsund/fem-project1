function poisson_3d_dirichlet(f, nr_of_mesh_nodes)
% The function takes a f function and a number of mesh nodes. It gives an
% approximation of u for the poisson in 3d problem.

% Get the triangulation elements, the vertices and edges
[p, tri, edge] = getSphere(nr_of_mesh_nodes);

% Get the stiffness matrix and the b vector before they are altered from
% the boundary conditions.
[A, b] = get_stiffness_matrix_and_load_vector_3D(nr_of_mesh_nodes, f, p, tri);

% Implementing Dirichlet boundary conditions
% Setting all columns and rows with an index on the edge to 0. Setting the
% diagonal elements to 1 for the zero rows. Let the elements in the
% b-vector be equal to 0 for the edge indices.

boundary = unique(edge);
A(boundary, :) = 0;
A(boundary, boundary) = eye(length(boundary));
b(boundary) = 0;

% Solving the system to find the weights for u.
u = A\b;

% Plot of the ball using triangulation and tetramesh
TR = triangulation(tri, p(:,1), p(:,2), p(:,3));
tetramesh(TR)
xlabel('x')
ylabel('y')
zlabel('z')

% u_analytical = zeros(nr_of_mesh_nodes, 1);
% k = @(x) sin(2 * pi * (x(1)^2 + x(2)^2 + x(3)^2));
% for i = 1 : nr_of_mesh_nodes
%     u_analytical(i) = k(p(i,:));
% end
% error = abs(u - u_analytical);
% max_error = max(error);
% 
% F = scatteredInterpolant(p(:,1), p(:,2), p(:,3), full(u));
% [X, Y, Z] = meshgrid(-1:0.01:0, -1:0.01:1, -1:0.01:1);
% R = (X.^2 + Y.^2 +Z.^2)<1;
% V = F(X,Y,Z).*R;

% % Plot of the numerical solution
% figure
% axis equal
% colormap parula
% colorbar
% camlight
% str = sprintf('Numerical solution using %d points. Max error: %f ', nr_of_mesh_nodes, max_error);
% title(str)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% for isovalue = min(u) : 0.05 : max(u)
%     isosurface(X, Y, Z, V, isovalue);
%     hold on
% end

% % Plot of analytical solution
% G = scatteredInterpolant(p(:,1), p(:,2), p(:,3), full(u_analytical));
% V_analytical = G(X,Y,Z).*R;
% 
% figure
% axis equal
% colormap parula
% colorbar
% camlight
% str = sprintf('Numerical solution using %d points.', nr_of_mesh_nodes);
% title(str)
% xlabel('x')
% ylabel('y')
% zlabel('z')
% for isovalue = min(u_analytical) : 0.05 : max(u_analytical)
%     isosurface(X, Y, Z, V_analytical, isovalue);
%     hold on
% end


end