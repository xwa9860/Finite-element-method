function [A, error] = CG_solver_KE_table (precond_KE_table, T, R, tol)

%{
find the a for K a=R
%A = K\R;  
%}
KE_table = precond_KE_table;
N = size(KE_table, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preconditioning acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%KE_table= precond_KE_table(KE_table)%T array Ne*1

%K = transpose(T) * K * T;


%R = transpose(T) * R;
for i=1:N
    R(i)= R(i)*T(i);
end
%display(R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CG solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial guess of A
A = zeros(N+1, 1);
r = -KE_multi(KE_table, A) + R;
Z = - KE_multi(KE_table, A) + R;
lambda = transpose(Z) * r / (transpose(Z) * KE_multi(KE_table, Z));
A = A + lambda * Z;

error = 1;
while error > tol 
    r = -KE_multi(KE_table, A) + R;
    theta = - transpose(r) * KE_multi(KE_table , Z)/ (transpose(Z)*KE_multi(KE_table, Z));
    Z = r + theta * Z;
    lambda = transpose(Z) * r / (transpose(Z) * KE_multi(KE_table, Z));
    A_new = A + lambda * Z;
    error = abs(lambda)* sqrt( transpose(Z) * KE_multi(KE_table, Z)) /sqrt(transpose(A) * KE_multi(KE_table, A));
    A=A_new;
end

for i =1:N+1
    A(i) = A(i)* T(i);
end



