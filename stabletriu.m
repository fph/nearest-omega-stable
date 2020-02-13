function T = stabletriu(M)

T = triu(M);
T = T - diag(max(real(diag(T)),0));