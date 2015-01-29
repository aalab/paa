function x = convexCircle(K)
% CONVEXCIRCLE(K) generates K equally spaced points on a circle on a simplex
%
% Copyright Sohan Seth sohan.seth@hiit.fi

r = 1; % Radius
theta = (0:(2/K):2) * pi; theta = theta(1:end-1); % Equispaced K angles

A = [[0, sqrt(2) * r]; [2 * (1 + sqrt(2)) / sqrt(3) * r, -r]; ...
    [-2 * (1 + sqrt(2)) / sqrt(3) * r, -r]]; % Vertices of triangle
A = A';

y = [r * cos(theta); r * sin(theta)]; % Samples

a1 = A(:,1) - A(:,3);
a2 = A(:,2) - A(:,3);
Z = - inv([a1, a2]) * bsxfun(@minus, A(:,3), y);
Z(3,:) = 1 - sum(Z); % Projection in the convex hull

B = [[1 0 0]; [0 1 0]; [0 0 1]];
x = B * Z; % Projection on simplex

%plot3(x(1,:), x(2,:), x(3,:), '.'), grid on