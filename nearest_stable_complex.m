function [S, xcost, x] = nearest_stable_complex(M, maxiter, timemax)

% Computes the nearest complex Hurwitz stable matrix to a matrix A (or at least a local minimum)
%
% [B, e, t, Q, U] = nearest_stable_complex(A, maxiter, timemax)
%
% B = arg min_{X Hurwitz stable}  ||X-A||_F
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
%M = randn(n);
problem.M = stiefelcomplexfactory(n, n, 1);
problem.cost = @(Q) norm(tril(Q'*M*Q,-1), 'fro')^2 + norm(max(real(diag(Q'*M*Q)),0))^2;
problem.grad = @(Q) Q*gradient_f(Q'*M*Q);   % the gradient on rotationmanifold and that on stiefelmanifold are represented differently, so here we need to multiply by Q and in the rotations example we didn't.

warning('off', 'manopt:getHessian:approx');

%checkgradient(problem);

options.tolgradnorm = 1e-12;
[x, xcost, info, options] = trustregions(problem, [], options);

S = x*stabletriu(x'*M*x)*x';