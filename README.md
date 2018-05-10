# ME3220Honors
1D element analysis for axial deflection of a beam

'''matlab
function M = mass(n,rho,A,L)
M = zeros(size(n+1));
C = rho*A*L/n/420;
matrix_m = 2/3*diag(ones(n+1,1))+1/6*diag(ones(n,1),1)+1/6*diag(ones(n,1),-1);
matrix_m(1,1) = 1/3; 
matrix_m(n+1,n+1) = 1/3;
M = C.*matrix_m;
M = M(2:n+1,2:n+1);
end
'''

'''matlab
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
end'''

'''matlab
function [w,u]=KM_solve(K,M)
[u,w] = eig(K,M);
w = sqrt(w);
w = diag(w);
end'''

'''matlab
A = 0.2;
E = 200*10^9;
rho = 8050;
L = 10;
for n = [4, 20, 50]
    K = stiffness(n,A,E,L);
    m = mass(n,rho,A,L);
    [w,u]=KM_solve(K,m);
    t_end = 2*pi/min(w);
    t = 0:0.00001:0.008;
    init = zeros(2*n,1);
    init(2*n) = -0.1;
    matrix_1 = zeros(n,2*n);
    matrix_2 = zeros(n,2*n);
    for i=2:n+1
        matrix_1(:,2*(i-1)) = u(:,i-1);
        matrix_2(:,2*i-3) = u(:,i-1);
    end
    matrix_3 = [matrix_1;matrix_2];
    coeff = matrix_3\init;
    solution = u(:,1)*coeff(1)*sin(w(1)*t);
    
    for i = 2:numel(w)
        solution = solution + u(:,i)*coeff(i)*sin(w(i)*t);
    end
    f = figure(1);
    plot(t*1000,solution(n,:)*1000)
    hold on;
    
end
title('Continuous Rod Deflection Over Time');
xlabel('Time (ms)')
ylabel('Position (mm)')
legend('n=4','n=20','n=50'); 
saveas(f,'HonorsProject.png');'''
