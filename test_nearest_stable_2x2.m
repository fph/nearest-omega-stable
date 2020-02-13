rng('default');
for trie = 1:100
    
    M = randn(2);
    
    [S, dist] = nearest_stable_2x2(M);
    [SM, distM] = manopt_nearest_stable_complex(M);
    if abs(dist - distM) > 1e-6
        keyboard
    end
end

% orthogonals are evil, so we test them separately

for trie = 1:100
    
    M = orth(randn(2));
    
    [S, dist] = nearest_stable_2x2(M);
    [SM, distM] = manopt_nearest_stable_complex(M);
    if abs(dist - distM) > 1e-6
        keyboard
    end
end
% ... with both signs for the determinant
for trie = 1:100
    
    M = orth(randn(2));
    M = diag([1,-1]) * M;
    
    [S, dist] = nearest_stable_2x2(M);
    [SM, distM] = manopt_nearest_stable_complex(M);
    if abs(dist - distM) > 1e-6
        keyboard
    end
end

