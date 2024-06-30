function [x,lambda] = eigSolver(A, B, epsilon, mu, x)
A = A-mu*B;% eigenvalue shift
dA = decomposition(A,'ldl','upper');%decomposition
x = x / norm(x);
u = B*x;
y = dA\u;
err = norm(y - x) / norm(y);
while err > epsilon
    x = y / norm(y);
    u = B*x;
    y = dA\u;
    lambda = x'*A*x/(x'*B*x);
    err = norm(y*lambda - x) ;
end
end