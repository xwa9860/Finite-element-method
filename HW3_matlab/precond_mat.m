function T = precond_mat(K)
%compute the diagnal matrix for CG preconditioning
N = size(K, 1);
t= zeros(N,1) ;
for i = 1 : N
    t(i) = 1/sqrt(K(i,i));
end
T = diag(t);


%{
test:

K = [1,2; 3,4]
T = precond_mat(K)
%T=[1, 0.5]


%}