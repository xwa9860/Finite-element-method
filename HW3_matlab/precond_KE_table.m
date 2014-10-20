function T = precond_KE_table(KE_table)
N=size(KE_table,1);
T = zeros(N+1, 1);

for i = 1: N
    T(i) = T(i)+ KE_table(i,1);
    T(i+1) = T(i+1) + KE_table(i,3);
end
