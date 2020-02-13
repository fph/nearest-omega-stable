function [S, e, t, x,U] = nearest_stable_real(M, maxiter, timemax)

% Computes the nearest real Hurwitz stable matrix to a real matrix A (or at least a local minimum)
%
% [B, e, t, Q, U] = nearest_stable_real(A, maxiter, timemax)
%
% B = arg min_{X Hurwitz stable}  ||X-A||_F
%
% B has quasi-Schur form B = Q*U*Q'.
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
   [U, L] = real_decomposition(Q'*M*Q);
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

[U, L] = real_decomposition(x'*M*x);
S = x*U*x';

infotable = struct2table(info);
e = sqrt(infotable.cost);
t = infotable.time;

end