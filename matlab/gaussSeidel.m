function [ x, err ] = gaussSeidel( x, A, b, iterations, omega )
%GAUSSSEIDEL Gauss-Seidel iteration
%   x - input initial point and output final point
% err - output error norm(A*x_i-b)
%   A - input matrix
%   b - input right-hand side
% iterations - how many iterations to do
% omega - Damping factor [0,2].

err = zeros([iterations,1]);

L = tril(A);
U = triu(A,1);

for i=1:iterations
    err(i) = norm(A*x-b);
    %x = L\(b-U*x);
    x = (1-omega)*x + omega*(L\(b-U*x));
    %fprintf(1,'%d %.2e\n',i,err(i));
end

end

