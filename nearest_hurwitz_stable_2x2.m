function [S, dist2] = nearest_hurwitz_stable_2x2(M)
% compute the nearest real Hurwitz stable matrix to a given 1x1 or 2x2 matrix

if size(M) == [1,1]
    S = M - max(0, real(M));
    dist2 = max(0, real(M))^2;
    return
end

assert(all(size(M)==[2,2]));

dist2 = inf;

% first attempt: the closest point on Tr<=0 (this may be equal to M if it is already stable)

S1 = M - 1/2*max(0,trace(M))*eye(2);
if det(S1) >= 0
    S = S1;
    dist2 = norm(M-S, 'fro')^2;
end

% second attempt: the closest point on Det=0

[U, Sigma, V] = svd(M);
S2 = U(:,1)*Sigma(1,1)*V(:,1)';
if trace(S2) <= 0
    if norm(M-S2, 'fro')^2 < dist2
       S = S2;
       dist2 = norm(M-S2, 'fro')^2;
    end
end

if dist2 == inf
% none worked, we must project on Tr=Det=0.
% We do this by first finding a basis s.t. M(1,1) == M(2,2)
    params = roots([M(2,2)-M(1,1) -2*(M(1,2)+M(2,1)) M(1,1)-M(2,2)]);
    x = params(end);
    Q=[1 x;-x 1]/sqrt(1+x^2);
    Z = Q'*M*Q;
    if abs(Z(2,1)) > abs(Z(1,2))
        S = Q(:,2)*Z(2,1)*Q(:,1)';
    else
        S = Q(:,1)*Z(1,2)*Q(:,2)';
    end
    dist2 = norm(M-S, 'fro')^2;
end
