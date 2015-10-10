function plot_three_meshes(mesh_sizes)

[p1, tri1] = getDisk(mesh_sizes(1));
[p2, tri2] = getDisk(mesh_sizes(2));
[p3, tri3] = getDisk(mesh_sizes(3));

subplot(1,3,1)
triplot(tri1, p1(:, 1), p1(:, 2))
str = sprintf('Disk %d points', mesh_sizes(1));
title(str)

subplot(1,3,2)
triplot(tri2, p2(:, 1), p2(:, 2))
str = sprintf('Disk %d points', mesh_sizes(2));
title(str)

subplot(1,3,3)
triplot(tri3, p3(:, 1), p3(:, 2))

str = sprintf('Disk %d points', mesh_sizes(3));
title(str)

end