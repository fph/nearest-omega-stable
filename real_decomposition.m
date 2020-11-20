function [U, L] = real_decomposition(M, f)
% computes a M = L+U decomposition with eigenvalues of U in a prescribed
% location
%
% f is a function that computes the projection for 1x1 or 2x2 blocks

if not(exist('f', 'var')) || (ischar(f) && strcmp(f, 'hurwitz'))
    f = @nearest_hurwitz_stable_2x2;
end

if ischar(f) && strcmp(f, 'schur')
    f = @nearest_schur_stable_2x2;
end
    
U = triu(M);
for i = 1:2:length(M)
    if i == length(M)
        I = i;
    else
        I = [i, i+1];
    end
    U(I,I) = f(M(I,I));
end
L = M - U;
