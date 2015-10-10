function [A, b] = poisson_2d_mixed_boundaries(f, nr_of_mesh_nodes)


% Get the triangulation elements, the vertices and edges
[p, tri, edge] = getDisk(nr_of_mesh_nodes);

% Get the stiffness matrix and the b vector before they are altered from
% the boundary conditions.
[A, b] = get_stiffness_matrix_and_b(nr_of_mesh_nodes, f, p, tri);

% Find all edge points, and then find sigma_N, i.e. the edge points at the
% top half of the disk.
all_edge_points = p(edge(:, 1),:);
sigma_N = edge(all_edge_points(:, 2) > 0, :);

g = @(x) 4 * pi;

% Alter the b-vector because of the new term from the Neumann boundary
% condition.
for k = 1:length(sigma_N)
    edge_point_1 = sigma_N(k, 1);
    edge_point_2 = sigma_N(k, 2);
    line_start = p(edge_point_1, :);
    line_end = p(edge_point_2, :);
    h = @(x) f(x) * g(x);
    quadLine2D(line_start, line_end, 4, h)
    b(edge_point_1) = b(edge_point_1) + quadLine2D(line_start, line_end, 4, h);
end


end