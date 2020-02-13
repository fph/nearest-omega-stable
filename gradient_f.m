function G = gradient_f(M)
% computes gradient of norm(ctril(Q'*M*Q), 'fro')^2 in Q = I

n = length(M);
L = tril(M, -1);
U = triu(M, 1);
L(1:n+1:end) = max(0,real(diag(M)));
U(1:n+1:end) = diag(M) - diag(L);

G = U*L' - L'*U;

G = G - G';