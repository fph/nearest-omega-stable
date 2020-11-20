function [S, e, t, x,U] = manopt_nearest_stable_complex_generic(M, fun, maxiter, timemax, x0)
% Compute the nearest (complex) Omega-stable matrix to a matrix A (or at least a local minimum)
%
% [B, e, t, Q, U] = manopt_nearest_stable_complex_generic(A, fun, maxiter, timemax, x0)
%
% B = arg min_{X stable}  ||X-A||_F
%
% B has quasi-Schur form B = Q*U*Q'.
% e, t contain the distance and elapsed time after each iteration (for diagnostic
% purposes).
% "fun" is either the string 'hurwitz', 'schur', or a function that computes the projection
% on the stable set Omega for complex values.
%
% Requires Manopt.

n = length(M);

problem.M = stiefelcomplexfactory(n, n, 1);

if not(exist('maxiter', 'var'))
    maxiter = 1000;
end
if not(exist('timemax', 'var'))
    timemax = 1000;
end
if not(exist('x0', 'var'))
    x0 = [];
end


function [f, store] = cost(Q, store)
    if not(isfield(store, 'L'))
        [store.U, store.L] = complex_decomposition(Q'*M*Q, fun);
    end
    f = norm(store.L, 'fro')^2;
end

function [g, store] = grad(Q, store)
    if not(isfield(store, 'L'))
        [store.U, store.L] = complex_decomposition(Q'*M*Q, fun);
    end
    if not(isfield(store, 'g'))
      store.g = store.U*store.L' - store.L'*store.U;
      store.g = store.g - store.g';
      store.g = Q*store.g;
    end
    g = store.g;
end

problem.cost = @cost;
problem.grad = @grad;
options.maxiter = maxiter;
options.maxtime = timemax;

warning('off', 'manopt:getHessian:approx');

%checkgradient(problem);
%return

options.tolgradnorm = 1e-8;

[x, xcost, info, options] = trustregions(problem, x0, options);

[U, L] = complex_decomposition(x'*M*x, fun);
S = x*U*x';

infotable = struct2table(info);
e = sqrt(infotable.cost);
t = infotable.time;

end