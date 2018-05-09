function M = mass(n,rho,A,L)
M = zeros(size(n+1));
C = rho*A*L/n/420;
matrix_m = 2/3*diag(ones(n+1,1))+1/6*diag(ones(n,1),1)+1/6*diag(ones(n,1),-1);
matrix_m(1,1) = 1/3; 
matrix_m(n+1,n+1) = 1/3;
M = C.*matrix_m;
M = M(2:n+1,2:n+1);
end



