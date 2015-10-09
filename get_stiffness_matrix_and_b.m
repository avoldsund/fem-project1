function [A, b] = get_stiffness_matrix_and_b(nr_of_mesh_nodes, f, vertices, triangles)

% Create an empty A-matrix and b-vector. These will be populated during the
% following for loop.
A = zeros(nr_of_mesh_nodes);
b = zeros(nr_of_mesh_nodes, 1);

% Populate the A-matrix and b-vector.
% Iterate over all elements (triangles). Find the vertices of the element, and calculate the
% basis functions. Every element has three vertices, and three basis functions. The basis functions are linear functions and of the form
% linear functions of the form: ax + by + c. We find a, b and c by solving
% a linear system.
for k = 1:length(triangles)
    indices = triangles(k,:);
    p1 = vertices(indices(1),:);
    p2 = vertices(indices(2),:);
    p3 = vertices(indices(3),:);
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
        % For the b-vector. The phi function is of the form ax + by + c.
        phi = @(x) C(1,i) * x(1) + C(2,i) * x(2) + C(3,i);
        h = @(x) f(x) * phi(x);
        b(indices(i)) = b(indices(i)) + quadrature2D(p1, p2, p3, 4, h);
        
        for j = 1:3
        % For the A-matrix
            g = @(x) dot(grad(:, i), grad(:, j));
            A(indices(i), indices(j)) = A(indices(i), indices(j)) + quadrature2D(p1, p2, p3, 4, g);
        end
    end
    
end

end