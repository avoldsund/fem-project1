function F = poisson_3d_dirichlet(f, nr_of_mesh_nodes)
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
A(edge(:,1), :) = 0;
A(edge(:,1), edge(:,1)) = eye(length(edge));
b(edge(:,1)) = 0;

% Solving the system to find the weights for u.
u = A\b;

% % Plot of the ball using triangulation and tetramesh
% TR = triangulation(tri, p(:,1), p(:,2), p(:,3));
% tetramesh(TR)
% xlabel('x')
% ylabel('y')
% zlabel('z')


F = scatteredInterpolant(p(:,1), p(:,2), p(:,3), full(u))
%patch(isosurface(F.Points(:,1), F.Points(:,2), F.Points(:,3), F.Values, 20))


[x,y,z,v] = flow;
p = patch(isosurface(x,y,z,v,-3));
isonormals(x,y,z,v,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
camlight 
lighting gouraud

end