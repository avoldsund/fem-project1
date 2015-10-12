function poisson_2d_mixed_boundaries(f, nr_of_mesh_nodes)
% The function takes a f function and a number of mesh nodes. It
% approximates the solution for u, and plots this.
% The plotted solution shows u on unit disc with the mixed boundary conditions.

% Get the triangulation elements, the vertices and edges
[p, tri, edge] = getDisk(nr_of_mesh_nodes);

% Get the stiffness matrix and the b vector before they are altered from
% the boundary conditions.
[A, b] = get_stiffness_matrix_and_load_vector_2D(nr_of_mesh_nodes, f, p, tri);

% Find all edge points, and then find sigma_N, i.e. the edge points at the
% top half of the disk.
all_edge_points = p(edge(:, 1),:);
sigma_N = edge(all_edge_points(:, 2) > 0, :);
sigma_D = setdiff(edge, sigma_N, 'rows');
g = @(x) 4 * pi;

% Alter the b-vector because of the new term from the Neumann boundary
% condition.
for k = 1:length(sigma_N)
    edge_point_1 = sigma_N(k, 1);
    edge_point_2 = sigma_N(k, 2);
    line_start = p(edge_point_1, :);
    line_end = p(edge_point_2, :);
    
    K = [line_start(1) line_start(2); line_end(1) line_end(2)];
    d_1 = [1;0]; d_2 = [0;1];
    scalars_1 = K\d_1;
    scalars_2 = K\d_2;
    
    phi_1 = @(x) scalars_1(1) * x(1) + scalars_1(2) * x(2);
    phi_2 = @(x) scalars_2(1) * x(1) + scalars_2(2) * x(2);
    h_1 = @(x) phi_1(x) * g(x);
    h_2 = @(x) phi_2(x) * g(x);
    
    b(edge_point_1) = b(edge_point_1) + quadLine2D(line_start, line_end, 4, h_1);
    b(edge_point_2) = b(edge_point_2) + quadLine2D(line_end, line_start, 4, h_2);
end

% Changing the A-matrix to correspond with boundary conditions
edge_indices_with_dirchlet = union(sigma_D(:,1), sigma_D(:,2));
A(edge_indices_with_dirchlet, :) = 0;
A(edge_indices_with_dirchlet, edge_indices_with_dirchlet) = eye(length(edge_indices_with_dirchlet));
b(edge_indices_with_dirchlet) = 0;
u = A\b;

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
figure
trimesh(tri, p(:,1), p(:,2), error)
colorbar
str = sprintf('Error plot. Max error: %f ', max_error);
title(str)
xlabel('x')
ylabel('y')

end