function [U, L] = complex_decomposition(M, f)
% computes a M = L+U decomposition with eigenvalues of U in a prescribed
% location
%
% f is the projection function on Omega

if not(exist('f', 'var')) || (ischar(f) && strcmp(f, 'hurwitz'))
    f = @(x) x - max(real(x), 0);
end

function y = schur_f(x)
    if abs(x) <= 1
        y = x;
    else
        y = x / abs(x);
    end
end

if ischar(f) && strcmp(f, 'schur')
    f = @schur_f;
end

n = length(M);
L = tril(M,-1);
U = triu(M, 1);
for i = 1:n
    U(i,i) = f(M(i,i));
    L(i,i) = M(i,i) - U(i,i);
end
end