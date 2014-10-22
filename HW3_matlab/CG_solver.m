function [A, error] = CG_solver (K, R, tol)

%{
find the a for K a=R
%A = K\R;  
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preconditioning acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T= precond_mat(K);
display(T);
K = transpose(T) * K * T;
display(R);
R = transpose(T) * R;
display(R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CG solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial guess of A
N = size(K, 1);
A = zeros(N, 1);
r = -K * A + R
Z = - K * A + R
lambda = transpose(Z) * r / (transpose(Z) * K * Z)
A = A + lambda * Z

error = 1;
while error > tol 
    r = -K*A + R
    theta = - transpose(r) * K * Z/ (transpose(Z)*K * Z)
    Z = r + theta * Z
    lambda = transpose(Z) * r / (transpose(Z) * K * Z)
    A_new = A + lambda * Z
    error = abs(lambda)* sqrt( transpose(Z) * K * Z) /sqrt(transpose(A) * K * A);
    A=A_new;
end
A = T*A;


