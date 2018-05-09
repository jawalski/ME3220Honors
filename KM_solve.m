function [w,u]=KM_solve(K,M)
[u,w] = eig(K,M);
w = sqrt(w);
w = diag(w);
end
