function K = stiffness(n,A,E,L)
K = zeros(size(n+1));
r_2 = A/pi();
I = (pi()/4)*r_2^2;
C = E*I/(L/n)^3;
%C = n*A*E/L;
matrix_k = 2*diag(ones(n+1,1))-1*diag(ones(n,1),1)-1*diag(ones(n,1),-1);
matrix_k(1,1) = 1;
matrix_k(n+1,n+1) = 1;
K = C.*matrix_k;
K = K(2:n+1,2:n+1);
end
