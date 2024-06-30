function [eigval,eigfun] = InvVec(A,B,lambda)
At = A-lambda*B;
eps = 1e-10;
eps1 = eps+1;
v0 = rand(size(A,1),1);
while eps1>eps
    x = At\(B*v0);
    eigval = 1/norm(x);
    v1 = x*eigval;
    eps1 = norm(v1-v0);
    v0 = v1;
end
eigval = eigval+lambda;
eigfun = v1;