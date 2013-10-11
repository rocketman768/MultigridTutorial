%% This script loads the  Laplace problem specified by:

% Size of variable x
N = 32;
% Normalized frequency [0,N) of initial guess x0
f = 1;

%% Do not edit below this line.
tmpB = zeros(N,3);

% Simple Laplacian
tmpB(:,1) = -ones(N,1);
tmpB(1:end,2) = 2*ones(N,1);
tmpB(:,3) = -ones(N,1);

% Without boundary effects
%tmpB(:,1) = -ones(N,1);
%tmpB(2:end-1,2) = 2*ones(N-2,1); tmpB(1,2) = 1; tmpB(end,2) = 1;
%tmpB(:,3) = -ones(N,1);

% Create sparse matrix
A = spdiags(tmpB,[-1,0,1],N,N);
% Create rhs
b = zeros(N,1);
% Create initial solution x0
t = linspace(0+1/N,1-1/N,N)';
x0 = sin( pi*f*t );

tmpn = N;
% NOTE: always doing the same number of small grids does NOT result in
% convergence that is independent of N! This only happens when the
% smallest grid is the same size, so the number of levels is O(log2(N)).
%up = cell(4,1);
%dowN = cell(4,1);
%for i=4:-1:1
% Always have the smallest grid be size 4.
up = cell(log2(N)-2,1);
dowN = cell(log2(N)-2,1);
for i=log2(N)-2:-1:1
    tmpijk = zeros(3*tmpn/2,3);
    % Beginning index.
    tmpijk(1:2,:) = [[1;1],[1;2],[1;0.5]];
    ijkndx = 3;
    bigndx = 3;
    for smallndx=2:tmpn/2
        tmpijk(ijkndx:ijkndx+2,:) = [repmat(smallndx,3,1),(bigndx-1:bigndx+1)',[0.5,1,0.5]'];
        ijkndx = ijkndx + 3;
        bigndx = bigndx + 2;
    end
    dowN{i} = sparse(tmpijk(1:ijkndx-1,1), tmpijk(1:ijkndx-1,2), tmpijk(1:ijkndx-1,3), tmpn/2, tmpn);
    up{i} = dowN{i}';
    tmpn = tmpn/2;
end

clear tmpB tmpn tmpijk ijkndx bigndx;