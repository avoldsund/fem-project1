function [A] = poisson_2d_dirichlet(f, nr_of_mesh_nodes)
%

[p, tri, edge] = getDisk(nr_of_mesh_nodes);
A = zeros(nr_of_mesh_nodes);


% Iterate over all elements, find the points of the element, and find the
% basis functions. The basis functions are linear functions and of the form
% ax + by + c. We find a, b and c in the foor loop.
for k = 1:length(tri)
    indices = tri(k,:);
    p1 = p(indices(1),:);
    p2 = p(indices(2),:);
    p3 = p(indices(3),:);
    K = [p1(1) p1(2) 1; p2(1) p2(2) 1; p3(1) p3(2) 1];
    b1 = [1;0;0];
    b2 = [0;1;0];
    b3 = [0;0;1];
    phi1_scalars = K\b1;
    phi2_scalars = K\b2;
    phi3_scalars = K\b3;
    C = [phi1_scalars, phi2_scalars, phi3_scalars];
    grad = C(1:2, :);
    for i = 1:3
        for j = 1:3
            g = @(x) dot(grad(:, i), grad(:, j));
            A(indices(i), indices(j)) = A(indices(i), indices(j)) + quadrature2D(p1, p2, p3, 4, g);
        end
    end
end


end