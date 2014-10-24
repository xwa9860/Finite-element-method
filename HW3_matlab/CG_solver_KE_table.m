function [A, error, it] = CG_solver_KE_table (KE_table, R, tol)

%{
find the a for K a=R
%A = K\R;  
%}

N = size(KE_table, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%preconditioning acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T= precond_KE_table(KE_table);%T array Ne*1

%K = transpose(T) * K * T;
for i=1:N
        KE_table(i,1) = T(i) * KE_table(i,1) *T(i); 
        KE_table(i,2) = T(i) * KE_table(i,2) *T(i+1); 
        KE_table(i,3) = T(i+1) * KE_table(i,3) *T(i+1); 
end

display(R);
%R = transpose(T) * R;
for i=1:N+1
    R(i)= R(i)*T(i);
end
display(R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CG solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial guess of A
A = zeros(N+1, 1);
r = -KE_multi(KE_table, A) + R;
Z = - KE_multi(KE_table, A) + R;
lambda = transpose(Z) * r / (transpose(Z) * KE_multi(KE_table, Z));
A = A + lambda * Z;
it=1;
error = 1;
while error > tol 
    r = -KE_multi(KE_table, A) + R;
    theta = - transpose(r) * KE_multi(KE_table , Z)/ (transpose(Z)*KE_multi(KE_table, Z));
    Z = r + theta * Z;
    lambda = transpose(Z) * r / (transpose(Z) * KE_multi(KE_table, Z));
    A_new = A + lambda * Z;
    error = abs(lambda)* sqrt( ( transpose(Z) * KE_multi(KE_table, Z)) /(transpose(A) * KE_multi(KE_table, A)));
    A=A_new;
    it=it+1;
end

for i =1:N+1
    A(i) = A(i)* T(i);
end



