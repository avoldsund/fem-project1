function [A, b] = get_stiffness_matrix_and_load_vector_3D(nr_of_mesh_nodes, f, vertices, tetrahedrons)

% Create an empty A-matrix and b-vector. These will be populated during the
% following for loop.
A = sparse(zeros(nr_of_mesh_nodes));
b = sparse(zeros(nr_of_mesh_nodes, 1));

for k = 1:length(tetrahedrons)
    indices = tetrahedrons(k,:);
    p1 = vertices(indices(1),:);
    p2 = vertices(indices(2),:);
    p3 = vertices(indices(3),:);
    p4 = vertices(indices(4),:);
    K = [p1(1) p1(2) p1(3) 1; p2(1) p2(2) p2(3) 1; p3(1) p3(2) p3(3) 1; p4(1) p4(2) p4(3) 1];
    b1 = [1;0;0;0];
    b2 = [0;1;0;0];
    b3 = [0;0;1;0];
    b4 = [0;0;0;1];
    phi1_scalars = K\b1;
    phi2_scalars = K\b2;
    phi3_scalars = K\b3;
    phi4_scalars = K\b4;
    C = [phi1_scalars, phi2_scalars, phi3_scalars, phi4_scalars];
    grad = C(1:3, :);
    
    for i = 1:4
        % For the b-vector. The phi function is of the form ax + by + cz + d.
        phi = @(x) C(1,i) * x(1) + C(2,i) * x(2) + C(3,i) * x(3) + C(4,i);
        h = @(x) f(x) * phi(x);
        b(indices(i)) = b(indices(i)) + quadrature3D(p1, p2, p3, p4, 5, h);
        
        for j = 1:4
        % For the A-matrix
            g = @(x) dot(grad(:, i), grad(:, j));
            A(indices(i), indices(j)) = A(indices(i), indices(j)) + quadrature3D(p1, p2, p3, p4, 5, g);
        end
    end
    
end

end