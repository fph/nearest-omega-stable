function [U, L] = real_decomposition(M, f);
% computes a M = L+U decomposition with eigenvalues of U in a prescribed
% location
%
% f is a function that 'splits' a 1x1 or 2x2 block into S + difference

if not(exist('f', 'var'))
    f = @nearest_stable_2x2;
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
