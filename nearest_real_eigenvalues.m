function [S, e, t, x,U] = nearest_real_eigenvalues(M, maxiter, timemax)

% Computes the nearest matrix with real eigenvalues to a real matrix A (or at least a local minimum)
%
% [B, e, t, Q, U] = nearest_real_eigenvalues(A, maxiter, timemax)
%
% B = arg min_{X has real eigenvalues}  ||X-A||_F
%
% B has Schur form B = Q*U*Q'.
% e, t contain the distance and time after each iteration (for diagnostic
% purposes).

if not(exist('maxiter', 'var'))
    maxiter = inf;
end
if not(exist('timemax', 'var'))
    timemax = inf;
end

n = length(M);
problem.M = rotationsfactory(n);

function [cost, grad] = costgrad(Q)
   A = Q'*M*Q;
   U = triu(A);
   L = tril(A, -1);
   cost = norm(L, 'fro')^2;
   grad = U*L' - L'*U;
   grad = grad - grad';
end

problem.costgrad = @costgrad;
options.maxiter = maxiter;
options.maxtime = timemax;

warning('off', 'manopt:getHessian:approx');

% checkgradient(problem);

options.tolgradnorm = 1e-12;

[x, xcost, info, options] = trustregions(problem, [], options);

U = triu(x'*M*x);
S = x*U*x';

infotable = struct2table(info);
e = sqrt(infotable.cost);
t = infotable.time;

end
