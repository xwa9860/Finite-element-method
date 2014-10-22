function [KE_table_cond, T] = precond_KE_table(KE_table)
N=size(KE_table,1);
T = zeros(N+1, 1);

for i = 1: N
    T(i) = T(i)+ (KE_table(i,1));
    T(i+1) = T(i+1) + (KE_table(i,3));
end

for i = 1:(N+1)
    T(i)=1/sqrt(T(i));
end

for i=1:N
        KE_table(i,1) = T(i) * KE_table(i,1) *T(i); 
        KE_table(i,2) = T(i) * KE_table(i,2) *T(i+1); 
        KE_table(i,3) = T(i+1) * KE_table(i,3) *T(i+1); 
end
KE_table_cond = KE_table;