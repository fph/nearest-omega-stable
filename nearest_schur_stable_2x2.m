function [B, dist2, type] = nearest_schur_stable_2x2(M)
% compute the nearest real Schur stable matrix to a given 1x1 or 2x2 matrix

% 1x1 case
if size(M) == [1,1]
    B = M;
    if abs(B) > 1
        B = B / abs(B);
    end
    dist2 = abs(M-B)^2;
    type = 1;
    return
end

assert(all(size(M)==[2,2]));

dist2 = inf;

% already stable
if all(abs(eig(M)) <= 1)
    B = M;
    dist2 = 0;
    type = 2;
    return
end

% B0
[U, S, V] = svd(M);
r = hyperbola_critical(diag(S));
for i = 1:length(r)
    B0 = U*diag([r(i); 1/r(i)])*V';
    if abs(trace(B0)) <= 1 + det(B0)
        if norm(B0 - M, 'fro')^2 < dist2
            B = B0;
            dist2 = norm(B0 - M, 'fro')^2;
            type = 3;
        end
    end
end

%Bplus
[U, S, V] = svd(M - eye(2));
Bplus = eye(2) + U(:,1)*S(1,1)*V(:,1)';
if det(Bplus) < 1 && -trace(Bplus) <= 1 + det(Bplus)
    if norm(Bplus - M, 'fro')^2 < dist2
        B = Bplus;
        dist2 = norm(Bplus - M, 'fro')^2;
        type = 4;
    end
end

%Bminus
[U, S, V] = svd(M + eye(2));
Bminus = -eye(2) + U(:,1)*S(1,1)*V(:,1)';
if det(Bminus) < 1 && trace(Bminus) <= 1 + det(Bminus)
    if norm(Bminus - M, 'fro')^2 < dist2
        B = Bminus;
        dist2 = norm(Bminus - M, 'fro')^2;
        type = 5;
    end
end

params = roots([M(2,2)-M(1,1) -2*(M(1,2)+M(2,1)) M(1,1)-M(2,2)]);
x = params(end);
G=[1 x;-x 1]/sqrt(1+x^2);
Z = G'*M*G;
% Z should have Z(1,1) == Z(2,2)

% B0plus/B0minus, always stable.
% We look directly for the closest matrix to Z
if trace(Z) >= 0
    B0pm = eye(2);
else
    B0pm = -eye(2);
end 
if abs(Z(1,2)) > abs(Z(2,1))
    B0pm(1,2) = Z(1,2);
else
    B0pm(2,1) = Z(2,1);
end
if norm(Z - B0pm, 'fro')^2 < dist2    
    B = G*B0pm*G';
    dist2 = norm(Z - B0pm, 'fro')^2;
    type = 6;
end

% Bpm
r = hyperbola_critical([Z(1,2), Z(2,1)]);
for i = 1:length(r)
    Bpm = [0 r(i); 1/r(i) 0];
    if norm(Z - Bpm, 'fro')^2 < dist2
        B = G*Bpm*G';
        dist2 = norm(Z - Bpm, 'fro')^2;
        type = 7;
    end
end
