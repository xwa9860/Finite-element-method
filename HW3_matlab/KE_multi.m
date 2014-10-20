function y = KE_multi(KE_table, q)
%q is a vector N*1
%calculate y = K * q
N = size(KE_table, 1);

y=zeros(N+1,1);

for e = 1:N
    Ke = KE_table(e, :);
    
    qe1 = q(e);
    qe2 = q(e+1);
    
    y(e) = y(e) + Ke(1)*qe1 + Ke(2)*qe2;
    y(e+1) = y(e+1) + Ke(2)*qe1 + Ke(3)*qe2; 
end